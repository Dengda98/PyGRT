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


def test_invalid_gf_source_and_freqband_and_dists():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(grn=HERE / "_tmp_args_grn", modelpath=MODEL)

        try:
            model.greenfn(depsrc=1.0, deprcv=0.0, dists=1.0, nt=8, dt=0.1, gf_source=["XX"])
        except ValueError as exc:
            assert "gf_source" in str(exc)
        else:
            raise AssertionError("invalid gf_source should raise ValueError")

        try:
            model.greenfn(depsrc=1.0, deprcv=0.0, dists=1.0, nt=8, dt=0.1, freqband=[1.0])
        except ValueError as exc:
            assert "freqband" in str(exc)
        else:
            raise AssertionError("short freqband should raise ValueError")

        try:
            model.greenfn(depsrc=1.0, deprcv=0.0, dists="10.5", nt=8, dt=0.1)
        except TypeError as exc:
            assert "dists" in str(exc)
        else:
            raise AssertionError("string dists should raise TypeError")

        try:
            model.greenfn(depsrc=1.0, deprcv=0.0, dists=[2.0, 1.0], nt=8, dt=0.1)
        except ValueError as exc:
            assert "strictly ascending" in str(exc)
        else:
            raise AssertionError("non-ascending dists should raise ValueError")
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_print_log_forwarded_to_run_grt():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(grn=HERE / "_tmp_args_grn", modelpath=MODEL)

        model.greenfn(depsrc=1.0, deprcv=0.0, dists=1.0, nt=8, dt=0.1, print_log=True)
        assert runner.kwargs[-1].get("print_log") is True
        assert "-s" not in runner.commands[-1]

        model.greenfn(depsrc=1.0, deprcv=0.0, dists=1.0, nt=8, dt=0.1, print_log=False)
        assert runner.kwargs[-1].get("print_log") is False
        assert "-s" in runner.commands[-1]
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_greenfn_default_and_optional_flags():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(grn=HERE / "_tmp_args_grn", modelpath=MODEL, topbound="free", botbound="halfspace")

        # 默认参数应显式带上 C 侧常用默认片段
        model.greenfn(depsrc=2.0, deprcv=0.0, dists=[1.0, 2.5], nt=32, dt=0.05)
        cmd = runner.commands[-1]
        assert_command_equals(
            cmd,
            [
                "greenfn",
                f"-M{MODEL}",
                "-D2/0",
                "-N32/0.05+w0.8+n1",
                "-R1,2.5",
                f"-O{model.grn}",
                "-BfH",
                "-H-1/-1",
                "-L0",
                "-K+k50+s2+e-1",
                "-E0/0",
            ],
        )

        # 各类可选参数拼接到正确的 CLI 选项
        model.greenfn(
            depsrc=3.5,
            deprcv=1.25,
            dists=10.0,
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
            f"-O{model.grn}",
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
        model.greenfn(
            depsrc=2.0,
            deprcv=0.0,
            dists=5.0,
            nt=16,
            dt=0.1,
            delayT0=0.5,
            ref_first_p=True,
            converg_method="DCM",
            safilonTol=1e-4,
            filonCut=1.0,
        )
        cmd = runner.commands[-1]
        assert_command_has(cmd, "-N16/0.1+w0.8+n1", "-L0+a0.0001+o1", "-Cd", "-Ep0.5")

        # 边界条件映射
        for top, bot, expected in [
            ("free", "halfspace", "-BfH"),
            ("rigid", "free", "-BrF"),
            ("halfspace", "rigid", "-BhR"),
        ]:
            model2 = pygrt.PyModel1D(grn=HERE / "_tmp_args_grn2", modelpath=MODEL, topbound=top, botbound=bot)
            model2.greenfn(depsrc=1.0, deprcv=0.0, dists=1.0, nt=8, dt=0.1)
            assert_command_has(runner.commands[-1], expected)

        # NONE 收敛方法
        model.greenfn(depsrc=1.0, deprcv=0.0, dists=1.0, nt=8, dt=0.1, converg_method="NONE")
        assert_command_has(runner.commands[-1], "-Cn")
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_static_greenfn_xy_and_dists():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(stgrn=HERE / "_tmp_args_static.nc", modelpath=MODEL)

        model.static_greenfn(
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
                "static_greenfn",
                f"-M{MODEL}",
                "-D2/3.3",
                f"-O{model.stgrn}",
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

        model.static_greenfn(depsrc=1.0, deprcv=0.0, dists=[0.0, 1.5, 3.0], safilonTol=1e-5, converg_method="PTAM")
        cmd = runner.commands[-1]
        assert_command_has(cmd, "static_greenfn", "-R0,1.5,3", "-L15+a1e-05", "-Cp", "-K+k50+e-1")
        assert "-X" not in " ".join(cmd)
        assert "-Y" not in " ".join(cmd)

        # 多深度：应拼出 -Ds/-Dr，且 stats 被忽略
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model.static_greenfn(depsrc=[1.0, 2.0, 3.0], deprcv=[0.0, 0.5], norths=[-2.0, 2.0, 1.0], easts=[-2.0, 2.0, 1.0], stats=True)
        assert any("stats" in str(w.message) for w in caught)
        cmd = runner.commands[-1]
        assert_command_has(cmd, "static_greenfn", "-Ds1,2,3", "-Dr0,0.5", "-X-2/2/1", "-Y-2/2/1")
        assert not any(str(tok) == "-S" for tok in cmd)
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_syn_source_and_time_function_options():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.pymod, runner)
    try:
        model = pygrt.PyModel1D(grn=HERE / "_tmp_args_grn", modelpath=MODEL)
        grn_root = HERE / "_tmp_args_grn"
        # 按 dist 匹配子目录，测试前需先准备假目录
        grn_dir = grn_root / f"{MODEL.name}_2_3.3_10"
        grn_dir.mkdir(parents=True, exist_ok=True)
        out = HERE / "_tmp_args_syn"

        # 爆炸源 + 时间函数 / 积分微分 / ZNE / upar
        model.syn(
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
                f"-G{grn_root}",
                "-A39.2",
                "-Su1e+24",
                f"-O{out / 'ex'}",
                "-R10",
                "-Dt/0.2/0.2/0.4",
                "-I1",
                "-J2",
                "-N",
                "-e",
            ],
        )

        # 单力源
        model.syn(dist=10.0, azimuth=12.0, scale=1e20, output_path=out / "sf", force=(2.0, -1.0, 4.0), time_function="t/0.1/0.3/0.6")
        assert_command_has(runner.commands[-1], "-F2/-1/4", "-Dt/0.1/0.3/0.6", "-S1e+20")

        # 剪切源 / 张裂源 / 矩张量
        model.syn(dist=10.0, azimuth=1.0, scale=1e22, output_path=out / "dc", strike=77.0, dip=88.0, rake=99.0, time_function="p/0.6")
        assert_command_has(runner.commands[-1], "-M77/88/99", "-Dp/0.6")

        model.syn(dist=10.0, azimuth=1.0, scale=1e22, output_path=out / "ts", strike=77.0, dip=88.0, time_function="p/0.6")
        assert_command_has(runner.commands[-1], "-M77/88", "-Dp/0.6")

        model.syn(
            dist=10.0,
            azimuth=1.0,
            scale=1e22,
            output_path=out / "mt",
            moment_tensor=(1.0, -2.0, -5.0, 0.5, 3.0, 1.2),
            time_function="r/3",
        )
        cmd = runner.commands[-1]
        assert_command_has(cmd, f"-G{grn_root}", "-T1/-2/-5/0.5/3/1.2", "-Dr/3")
    finally:
        _restore_run_grt(pygrt.pymod, original)


def test_source_type_is_inferred_from_source_parameters():
    import inspect

    assert "source" not in inspect.signature(pygrt.PyModel1D.syn).parameters
    assert "source" not in inspect.signature(pygrt.PyModel1D.static_syn).parameters
    assert "source" not in inspect.signature(pygrt.utils.okada).parameters

    runner = CapturedRunner()
    original_pymod = _patch_run_grt(pygrt.pymod, runner)
    original_utils = _patch_run_grt(pygrt.utils, runner)
    try:
        model = pygrt.PyModel1D(grn=HERE / "_tmp_args_grn", stgrn=HERE / "_tmp_args_stgrn.nc", modelpath=MODEL)
        model.syn(dist=10.0, azimuth=1.0, scale=1e20, output_path=HERE / "_tmp_args_inferred_sf", force=(1.0, 2.0, 3.0))
        assert_command_has(runner.commands[-1], "-F1/2/3")

        model.static_syn(scale=1e20, output_path=HERE / "_tmp_args_inferred_dc.nc", strike=33.0, dip=50.0, rake=120.0)
        assert_command_has(runner.commands[-1], "-M33/50/120")

        pygrt.utils.okada(
            modelparams=(6.0, 3.464, 2.7),
            depsrc=2.0,
            deprcv=0.0,
            norths=(-1.0, 1.0, 1.0),
            easts=(-1.0, 1.0, 1.0),
            output_path=HERE / "_tmp_args_okada.nc",
            scale=1e20,
            strike=33.0,
            dip=50.0,
        )
        assert_command_has(runner.commands[-1], "-M33/50")

        try:
            model.syn(dist=10.0, azimuth=1.0, scale=1e20, output_path=HERE / "_tmp_args_invalid", strike=33.0)
        except ValueError as exc:
            assert "strike and dip" in str(exc)
        else:
            raise AssertionError("incomplete source geometry should raise ValueError")

        try:
            model.syn(
                dist=10.0,
                azimuth=1.0,
                scale=1e20,
                output_path=HERE / "_tmp_args_mixed_source",
                force=(1.0, 2.0, 3.0),
                strike=33.0,
                dip=50.0,
            )
        except ValueError as exc:
            assert "mutually exclusive" in str(exc)
        else:
            raise AssertionError("mixed source parameters should raise ValueError")

        try:
            model.syn(dist=10.0, azimuth=1.0, scale=1e20, output_path=HERE / "_tmp_args_removed_source", source="EX")
        except TypeError as exc:
            assert "source" in str(exc)
        else:
            raise AssertionError("source should no longer be a public synthesis argument")
    finally:
        _restore_run_grt(pygrt.pymod, original_pymod)
        _restore_run_grt(pygrt.utils, original_utils)


def test_static_syn_and_tensor_postprocess_args():
    runner = CapturedRunner()
    original_pymod = _patch_run_grt(pygrt.pymod, runner)
    original_utils = _patch_run_grt(pygrt.utils, runner)
    try:
        model = pygrt.PyModel1D(stgrn=HERE / "_tmp_args_static.nc", modelpath=MODEL)
        out = HERE / "_tmp_args_static_syn.nc"

        model.static_syn(
            scale=1e24,
            output_path=out,
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
                "static_syn",
                f"-G{model.stgrn}",
                f"-O{out}",
                "-Su1e+24",
                "-M33/50/120",
                "-X-5/5/1",
                "-Y-4/4/2",
                "-N",
                "-e",
            ],
        )

        # 多深度点源 + 任意接收点
        model.static_syn(scale=1e20, output_path=out, depsrc=2.0, recv_points=HERE / "rcv.txt")
        assert_command_has(runner.commands[-1], "static_syn", f"-G{model.stgrn}", f"-O{out}", "-S1e+20", "-Ds2", f"-Q{HERE / 'rcv.txt'}")

        # 多深度点源 + 新网格 + 台站深度
        model.static_syn(scale=1e20, output_path=out, depsrc=2.0, deprcv=0.5, norths=[-2.0, 2.0, 1.0], easts=[-2.0, 2.0, 1.0])
        assert_command_has(runner.commands[-1], "static_syn", "-Ds2", "-Dr0.5", "-X-2/2/1", "-Y-2/2/1")

        # 有限断层
        model.static_syn(output_path=out, src_fault=HERE / "cfaults.inp", src_fault_size=(1.0, 2.0), calc_upar=True)
        assert_command_has(runner.commands[-1], "static_syn", f"-G{model.stgrn}", f"-O{out}", f"-C{HERE / 'cfaults.inp'}+i1/2", "-e")

        # 张量后处理：动态模块处理 SAC 目录，静态模块处理 NetCDF 文件
        dyn = HERE / "_tmp_tensor_dyn"
        dyn.mkdir(exist_ok=True)
        stc = HERE / "_tmp_tensor_static.nc"
        stc.write_bytes(b"placeholder")

        pygrt.utils.strain(dyn)
        assert_command_equals(runner.commands[-1], ["strain", str(dyn)])
        pygrt.utils.rotation(dyn)
        assert_command_equals(runner.commands[-1], ["rotation", str(dyn)])
        pygrt.utils.static_strain(stc)
        assert_command_equals(runner.commands[-1], ["static_strain", str(stc)])
        pygrt.utils.static_rotation(stc)
        assert_command_equals(runner.commands[-1], ["static_rotation", str(stc)])
        pygrt.utils.static_stress(stc)
        assert_command_equals(runner.commands[-1], ["static_stress", str(stc)])

        try:
            pygrt.utils.strain(stc)
        except ValueError:
            pass
        else:
            raise AssertionError("dynamic tensor modules should reject static files")

        try:
            pygrt.utils.static_strain(dyn)
        except ValueError:
            pass
        else:
            raise AssertionError("static tensor modules should reject dynamic directories")
    finally:
        _restore_run_grt(pygrt.pymod, original_pymod)
        _restore_run_grt(pygrt.utils, original_utils)


def test_static_sproj_and_coulomb_args():
    runner = CapturedRunner()
    original = _patch_run_grt(pygrt.utils, runner)
    static = HERE / "_tmp_args_sproj.nc"
    receiver = HERE / "_tmp_args_sproj_q.txt"
    dynamic = HERE / "_tmp_args_sproj_dir"
    try:
        static.write_bytes(b"placeholder")
        receiver.write_text("0 0 0 33 44 55\n", encoding="utf-8")
        dynamic.mkdir(exist_ok=True)

        pygrt.utils.static_sproj(static, strike=33.0, dip=44.0, rake=55.0)
        assert_command_equals(runner.commands[-1], ["static_sproj", f"-G{static}", "-M33/44/55"])

        pygrt.utils.static_sproj(static, rake=55.0, force_rake=True)
        assert_command_equals(runner.commands[-1], ["static_sproj", f"-G{static}", "-M55+f"])

        pygrt.utils.static_sproj(static, recv_points=receiver)
        assert_command_equals(runner.commands[-1], ["static_sproj", f"-G{static}", f"-Q{receiver}"])

        pygrt.utils.static_coulomb(static, 0.6)
        assert_command_equals(runner.commands[-1], ["static_coulomb", f"-G{static}", "-F0.6"])

        try:
            pygrt.utils.static_sproj(dynamic)
        except ValueError:
            pass
        else:
            raise AssertionError("dynamic synthesis directories should be rejected")
    finally:
        _restore_run_grt(pygrt.utils, original)
        static.unlink(missing_ok=True)
        receiver.unlink(missing_ok=True)
        dynamic.rmdir()


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

        st_strain = pygrt.utils.strain(dyn, return_result=True)
        assert_command_equals(runner.commands[-1], ["strain", str(dyn)])
        assert len(st_strain) == 1
        assert np.allclose(st_strain[0].data, samples["strain_EE.sac"])

        st_rot = pygrt.utils.rotation(dyn, return_result=True)
        assert_command_equals(runner.commands[-1], ["rotation", str(dyn)])
        assert len(st_rot) == 1
        assert np.allclose(st_rot[0].data, samples["rotation_NE.sac"])

        st_stress = pygrt.utils.stress(dyn, return_result=True)
        assert_command_equals(runner.commands[-1], ["stress", str(dyn)])
        assert len(st_stress) == 1
        assert np.allclose(st_stress[0].data, samples["stress_EE.sac"])
    finally:
        _restore_run_grt(pygrt.utils, original)
        shutil.rmtree(dyn, ignore_errors=True)


def test_renamed_interfaces_fail_with_migration_message():
    import inspect

    model = pygrt.PyModel1D()
    method_replacements = {
        "compute_travt1d": "travt",
        "compute_grn": "greenfn",
        "compute_static_grn": "static_greenfn",
        "compute_syn": "syn",
        "compute_static_syn": "static_syn",
    }
    for old_name, new_name in method_replacements.items():
        signature = inspect.signature(getattr(pygrt.PyModel1D, old_name))
        assert list(signature.parameters) == ["self", "args", "kwargs"]
        try:
            getattr(model, old_name)()
        except RuntimeError as exc:
            assert new_name in str(exc)
        else:
            raise AssertionError(f"{old_name} should reject the old interface")

    function_replacements = {
        "compute_okada": "okada",
        "compute_strain": "strain",
        "compute_rotation": "rotation",
        "compute_stress": "stress",
        "compute_sproj": "static_sproj",
        "compute_coulomb": "static_coulomb",
        "solve_lamb1": "lamb1",
    }
    for old_name, new_name in function_replacements.items():
        signature = inspect.signature(getattr(pygrt.utils, old_name))
        assert list(signature.parameters) == ["args", "kwargs"]
        try:
            getattr(pygrt.utils, old_name)()
        except RuntimeError as exc:
            assert new_name in str(exc)
        else:
            raise AssertionError(f"{old_name} should reject the old interface")


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
        test_invalid_gf_source_and_freqband_and_dists,
        test_print_log_forwarded_to_run_grt,
        test_greenfn_default_and_optional_flags,
        test_static_greenfn_xy_and_dists,
        test_syn_source_and_time_function_options,
        test_source_type_is_inferred_from_source_parameters,
        test_static_syn_and_tensor_postprocess_args,
        test_static_sproj_and_coulomb_args,
        test_tensor_return_result_reads_prefix_only,
        test_renamed_interfaces_fail_with_migration_message,
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
