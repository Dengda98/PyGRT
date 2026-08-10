"""
验证 Python API 到 grt CLI 参数的映射是否正确

不实际执行数值积分，只捕获 run_grt 收到的命令行参数
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pygrt
from compare_func import assert_command_equals, assert_command_has


HERE = Path(__file__).resolve().parent
MODEL = (HERE.parent / "milrow").resolve()


class CapturedRunner:
    """临时替换 run_grt，记录调用参数"""

    def __init__(self):
        self.commands = []
        self.kwargs = []

    def __call__(self, command, **kwargs):
        self.commands.append([str(item) for item in command])
        self.kwargs.append(kwargs)


def _patch_run_grt(monkey_target, runner: CapturedRunner):
    original = monkey_target.run_grt
    monkey_target.run_grt = runner
    return original


def _restore_run_grt(monkey_target, original):
    monkey_target.run_grt = original


def test_invalid_gf_source_and_freqband_and_distarr():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(MODEL)
        model.set_dynamic_grn_path(HERE / "_tmp_args_grn")

        try:
            model.compute_grn(depsrc=1.0, deprcv=0.0, distarr=1.0, nt=8, dt=0.1, gf_source=["XX"])
        except ValueError as exc:
            assert "gf_source" in str(exc)
        else:
            raise AssertionError("invalid gf_source should raise ValueError")

        try:
            model.compute_grn(depsrc=1.0, deprcv=0.0, distarr=1.0, nt=8, dt=0.1, freqband=[1.0])
        except ValueError as exc:
            assert "freqband" in str(exc)
        else:
            raise AssertionError("short freqband should raise ValueError")

        try:
            model.compute_grn(depsrc=1.0, deprcv=0.0, distarr="10.5", nt=8, dt=0.1)
        except TypeError as exc:
            assert "distarr" in str(exc)
        else:
            raise AssertionError("string distarr should raise TypeError")
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_print_log_forwarded_to_run_grt():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(MODEL)
        model.set_dynamic_grn_path(HERE / "_tmp_args_grn")

        model.compute_grn(depsrc=1.0, deprcv=0.0, distarr=1.0, nt=8, dt=0.1, print_log=True)
        assert runner.kwargs[-1].get("print_log") is True
        assert "-s" not in runner.commands[-1]

        model.compute_grn(depsrc=1.0, deprcv=0.0, distarr=1.0, nt=8, dt=0.1, print_log=False)
        assert runner.kwargs[-1].get("print_log") is False
        assert "-s" in runner.commands[-1]
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_compute_grn_default_and_optional_flags():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(MODEL, "free", "halfspace")
        model.set_dynamic_grn_path(HERE / "_tmp_args_grn")

        # 默认参数应显式带上 C 侧常用默认片段
        model.compute_grn(
            depsrc=2.0,
            deprcv=0.0,
            distarr=[1.0, 2.5],
            nt=32,
            dt=0.05,
        )
        cmd = runner.commands[-1]
        assert_command_equals(
            cmd,
            [
                "greenfn",
                f"-M{MODEL}",
                "-D2/0",
                "-N32/0.05+w0.8+n1",
                "-R1,2.5",
                f"-O{model.dynamic_grn_path}",
                "-BfH",
                "-H-1/-1",
                "-L0",
                "-K+k50+s2+e-1",
                "-E0/0",
            ],
        )

        # 各类可选参数拼接到正确的 CLI 选项
        model.compute_grn(
            depsrc=3.5,
            deprcv=1.25,
            distarr=10.0,
            nt=64,
            dt=0.02,
            upsampling_n=2,
            freqband=(0.1, 5.0),
            zeta=0.6,
            keepAllFreq=True,
            vmin_ref=1.5,
            keps=1e-3,
            ampk=3.0,
            k0=40.0,
            use_kmax_ref=True,
            Length=20.0,
            filonLength=5.0,
            filonCut=2.0,
            converg_method="PTAM",
            delayT0=1.2,
            delayV0=3.4,
            skipImagComps=True,
            calc_upar=True,
            gf_source=["EX", "DC", "HF"],
            statsidxs=[0, 3, 7],
            print_log=False,
        )
        cmd = runner.commands[-1]
        assert_command_has(
            cmd,
            "greenfn",
            f"-M{MODEL}",
            "-D3.5/1.25",
            "-N64/0.02+w0.6+n2+a+f",
            "-R10",
            f"-O{model.dynamic_grn_path}",
            "-BfH",
            "-H0.1/5",
            "-L20+l5+o2",
            "-Cp",
            "-K+k40+f+s3+e0.001+v1.5",
            "-E1.2/3.4",
            "-Gesh",
            "-S0,3,7",
            "-e",
            "-s",
        )

        # ref_first_p 对应 -Ep
        model.compute_grn(
            depsrc=2.0,
            deprcv=0.0,
            distarr=5.0,
            nt=16,
            dt=0.1,
            delayT0=0.5,
            ref_first_p=True,
            converg_method="DCM",
            safilonTol=1e-4,
            filonCut=1.0,
        )
        cmd = runner.commands[-1]
        assert_command_has(
            cmd,
            "-N16/0.1+w0.8+n1",
            "-L0+a0.0001+o1",
            "-Cd",
            "-Ep0.5",
        )

        # 边界条件映射
        for top, bot, expected in [
            ("free", "halfspace", "-BfH"),
            ("rigid", "free", "-BrF"),
            ("halfspace", "rigid", "-BhR"),
        ]:
            model2 = pygrt.PyModel1D(MODEL, top, bot)
            model2.set_dynamic_grn_path(HERE / "_tmp_args_grn2")
            model2.compute_grn(
                depsrc=1.0,
                deprcv=0.0,
                distarr=1.0,
                nt=8,
                dt=0.1,
            )
            assert_command_has(runner.commands[-1], expected)

        # NONE 收敛方法
        model.compute_grn(
            depsrc=1.0,
            deprcv=0.0,
            distarr=1.0,
            nt=8,
            dt=0.1,
            converg_method="NONE",
        )
        assert_command_has(runner.commands[-1], "-Cn")
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_compute_static_grn_xy_and_distarr():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(MODEL)
        model.set_static_grn_path(HERE / "_tmp_args_static.nc")

        model.compute_static_grn(
            depsrc=2.0,
            deprcv=3.3,
            norths=[-3.1, 3.1, 0.6],
            easts=[-4.1, 4.1, 0.8],
            calc_upar=True,
            stats=True,
            use_kmax_ref=True,
            k0=30.0,
            keps=1e-4,
            Length=12.0,
            filonLength=2.0,
            filonCut=0.5,
            converg_method="DCM",
        )
        cmd = runner.commands[-1]
        assert_command_equals(
            cmd,
            [
                "static",
                "greenfn",
                f"-M{MODEL}",
                "-D2/3.3",
                f"-O{model.static_grn_path}",
                "-BfH",
                "-X-3.1/3.1/0.6",
                "-Y-4.1/4.1/0.8",
                "-L12+l2+o0.5",
                "-Cd",
                "-K+k30+f+e0.0001",
                "-S",
                "-e",
            ],
        )

        model.compute_static_grn(
            depsrc=1.0,
            deprcv=0.0,
            distarr=[0.0, 1.5, 3.0],
            safilonTol=1e-5,
            converg_method="PTAM",
        )
        cmd = runner.commands[-1]
        assert_command_has(
            cmd,
            "static",
            "greenfn",
            "-R0,1.5,3",
            "-L15+a1e-05",
            "-Cp",
            "-K+k50+e-1",
        )
        assert "-X" not in " ".join(cmd)
        assert "-Y" not in " ".join(cmd)

        # 多深度：应拼出 -Ds/-Dr，且 stats 被忽略
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model.compute_static_grn(
                depsrc=[1.0, 2.0, 3.0],
                deprcv=[0.0, 0.5],
                norths=[-2.0, 2.0, 1.0],
                easts=[-2.0, 2.0, 1.0],
                stats=True,
            )
        assert any("stats" in str(w.message) for w in caught)
        cmd = runner.commands[-1]
        assert_command_has(cmd, "static", "greenfn", "-Ds1,2,3", "-Dr0,0.5", "-X-2/2/1", "-Y-2/2/1")
        assert not any(str(tok) == "-S" for tok in cmd)
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_compute_syn_source_and_time_function_options():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(MODEL)
        grn_root = HERE / "_tmp_args_grn"
        model.set_dynamic_grn_path(grn_root)
        # 按 dist 匹配子目录，测试前需先准备假目录
        grn_dir = grn_root / f"{MODEL.name}_2_3.3_10"
        grn_dir.mkdir(parents=True, exist_ok=True)
        out = HERE / "_tmp_args_syn"

        # 爆炸源 + 时间函数 / 积分微分 / ZNE / upar
        model.compute_syn(
            dist=10.0,
            azimuth=39.2,
            scale=1e24,
            output_path=out / "ex",
            time_function="t/0.2/0.2/0.4",
            integrate_order=1,
            differentiate_order=2,
            scale_with_mu=True,
            zne=True,
            calc_upar=True,
        )
        cmd = runner.commands[-1]
        assert_command_equals(
            cmd,
            [
                "syn",
                f"-G{grn_dir}",
                "-A39.2",
                "-Su1e+24",
                f"-O{out / 'ex'}",
                "-Dt/0.2/0.2/0.4",
                "-I1",
                "-J2",
                "-N",
                "-e",
            ],
        )

        # 单力源
        model.compute_syn(
            dist=10.0,
            azimuth=12.0,
            scale=1e20,
            output_path=out / "sf",
            source="SF",
            force=(2.0, -1.0, 4.0),
            time_function="t/0.1/0.3/0.6",
        )
        assert_command_has(
            runner.commands[-1],
            "-F2/-1/4",
            "-Dt/0.1/0.3/0.6",
            "-S1e+20",
        )

        # 剪切源 / 张裂源 / 矩张量
        model.compute_syn(
            dist=10.0,
            azimuth=1.0,
            scale=1e22,
            output_path=out / "dc",
            source="DC",
            strike=77.0,
            dip=88.0,
            rake=99.0,
            time_function="p/0.6",
        )
        assert_command_has(runner.commands[-1], "-M77/88/99", "-Dp/0.6")

        model.compute_syn(
            dist=10.0,
            azimuth=1.0,
            scale=1e22,
            output_path=out / "ts",
            source="TS",
            strike=77.0,
            dip=88.0,
            time_function="p/0.6",
        )
        assert_command_has(runner.commands[-1], "-M77/88", "-Dp/0.6")

        model.compute_syn(
            dist=10.0,
            azimuth=1.0,
            scale=1e22,
            output_path=out / "mt",
            source="MT",
            moment_tensor=(1.0, -2.0, -5.0, 0.5, 3.0, 1.2),
            time_function="r/3",
        )
        cmd = runner.commands[-1]
        assert_command_has(
            cmd,
            f"-G{grn_dir}",
            "-T1/-2/-5/0.5/3/1.2",
            "-Dr/3",
        )
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_compute_static_syn_and_tensor_postprocess_args():
    runner = CapturedRunner()
    original_pymod = _patch_run_grt(pygrt.pymod, runner)
    original_utils = _patch_run_grt(pygrt.utils, runner)
    try:
        model = pygrt.PyModel1D(MODEL)
        model.set_static_grn_path(HERE / "_tmp_args_static.nc")
        out = HERE / "_tmp_args_static_syn.nc"

        model.compute_static_syn(
            scale=1e24,
            output_path=out,
            source="DC",
            strike=33.0,
            dip=50.0,
            rake=120.0,
            norths=[-5.0, 5.0, 1.0],
            easts=[-4.0, 4.0, 2.0],
            zne=True,
            calc_upar=True,
            scale_with_mu=True,
        )
        assert_command_equals(
            runner.commands[-1],
            [
                "static",
                "syn",
                f"-G{model.static_grn_path}",
                "-Su1e+24",
                f"-O{out}",
                "-M33/50/120",
                "-X-5/5/1",
                "-Y-4/4/2",
                "-N",
                "-e",
            ],
        )

        # 张量后处理：目录走动态模块，文件走 static 模块
        dyn = HERE / "_tmp_tensor_dyn"
        dyn.mkdir(exist_ok=True)
        stc = HERE / "_tmp_tensor_static.nc"
        stc.write_bytes(b"placeholder")

        pygrt.utils.compute_strain(dyn)
        assert_command_equals(runner.commands[-1], ["strain", str(dyn)])
        pygrt.utils.compute_rotation(dyn)
        assert_command_equals(runner.commands[-1], ["rotation", str(dyn)])
        pygrt.utils.compute_stress(stc)
        assert_command_equals(
            runner.commands[-1],
            ["static", "stress", str(stc)],
        )
    finally:
        _restore_run_grt(pygrt.pymod, original_pymod)
        _restore_run_grt(pygrt.utils, original_utils)


def test_tensor_return_result_reads_prefix_only():
    """return_result=True 时只读对应前缀的 SAC，不与位移/其它张量混淆"""
    import shutil

    import numpy as np
    from obspy import Trace

    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.utils, runner)
    dyn = HERE / "_tmp_tensor_return"
    try:
        if dyn.exists():
            shutil.rmtree(dyn)
        dyn.mkdir()

        # 同 channel 名 EE/NE，只能靠文件名前缀区分
        samples = {
            "Z.sac": np.array([1.0, 2.0], dtype=np.float32),
            "strain_EE.sac": np.array([10.0, 20.0], dtype=np.float32),
            "stress_EE.sac": np.array([100.0, 200.0], dtype=np.float32),
            "rotation_NE.sac": np.array([3.0, 4.0], dtype=np.float32),
        }
        for name, data in samples.items():
            tr = Trace(data=data.copy())
            tr.stats.delta = 0.1
            tr.stats.channel = name.split("_")[-1].removesuffix(".sac") if "_" in name else name[0]
            tr.write(str(dyn / name), format="SAC")

        st_strain = pygrt.utils.compute_strain(dyn, return_result=True)
        assert_command_equals(runner.commands[-1], ["strain", str(dyn)])
        assert len(st_strain) == 1
        assert np.allclose(st_strain[0].data, samples["strain_EE.sac"])

        st_rot = pygrt.utils.compute_rotation(dyn, return_result=True)
        assert_command_equals(runner.commands[-1], ["rotation", str(dyn)])
        assert len(st_rot) == 1
        assert np.allclose(st_rot[0].data, samples["rotation_NE.sac"])

        st_stress = pygrt.utils.compute_stress(dyn, return_result=True)
        assert_command_equals(runner.commands[-1], ["stress", str(dyn)])
        assert len(st_stress) == 1
        assert np.allclose(st_stress[0].data, samples["stress_EE.sac"])
    finally:
        _restore_run_grt(pygrt.utils, original)
        shutil.rmtree(dyn, ignore_errors=True)


def test_format_helpers():
    from pygrt.cli import format_float, format_range

    assert format_float(1.0) == "1"
    assert format_float(1.25) == "1.25"
    assert format_float(1e24) == "1e+24"
    assert format_range([-3.1, 3.1, 0.6], "norths") == "-3.1/3.1/0.6"
    try:
        format_range([1, 2], "norths")
    except ValueError:
        pass
    else:
        raise AssertionError("format_range should reject non-3-length input")


def main():
    tests = [
        test_format_helpers,
        test_invalid_gf_source_and_freqband_and_distarr,
        test_print_log_forwarded_to_run_grt,
        test_compute_grn_default_and_optional_flags,
        test_compute_static_grn_xy_and_distarr,
        test_compute_syn_source_and_time_function_options,
        test_compute_static_syn_and_tensor_postprocess_args,
        test_tensor_return_result_reads_prefix_only,
    ]
    for func in tests:
        print(f"[RUN] {func.__name__}")
        func()
        print(f"[OK ] {func.__name__}")
    print("All CLI argument tests passed.")


if __name__ == "__main__":
    # 保证可直接脚本运行
    sys.path.insert(0, str(HERE))
    main()
