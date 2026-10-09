from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from obspy import read

import pygrt


MODEL_PARAMS = (8.0, 4.62, 3.3)
COMMON_PARAMS = {
    "modelparams": MODEL_PARAMS,
    "dist": 10.0,
    "nt": 16,
    "dt": 0.01,
    "azimuth": 30.0,
    "scale": 1e20,
    "print_log": False,
}


with TemporaryDirectory(prefix="pygrt-lamb-") as temporary:
    output_root = Path(temporary)

    try:
        pygrt.utils.lamb(
            **COMMON_PARAMS,
            depsrc=5.0,
            deprcv=0.0,
            output_path=output_root / "empty_phases",
            force=(0.5, -1.0, 2.0),
            phases=[],
        )
    except ValueError:
        pass
    else:
        raise ValueError("lamb should reject an empty phase list.")

    pygrt.utils.lamb(
        **COMMON_PARAMS,
        depsrc=0.0,
        deprcv=0.0,
        output_path=output_root / "surface",
        force=(0.5, -1.0, 2.0),
    )
    pygrt.utils.lamb(
        **COMMON_PARAMS,
        depsrc=5.0,
        deprcv=0.0,
        output_path=output_root / "lamb2",
        moment_tensor=(1.0, -2.0, 3.0, 0.5, 1.2, -0.7),
    )
    pygrt.utils.lamb(
        **COMMON_PARAMS,
        depsrc=5.0,
        deprcv=0.0,
        output_path=output_root / "lamb2_phases",
        moment_tensor=(1.0, -2.0, 3.0, 0.5, 1.2, -0.7),
        phases=["P", "S", "SP"],
    )
    pygrt.utils.lamb(
        **COMMON_PARAMS,
        depsrc=5.0,
        deprcv=1.0,
        output_path=output_root / "lamb3",
        strike=100.0,
        dip=30.0,
        rake=70.0,
    )

    # 西侧台站的 atan2 方位角为负，点文件和极坐标必须给出相同的非零波形
    points = output_root / "western_receiver.txt"
    points.write_text("0 -10 1\n")
    common = dict(modelparams=MODEL_PARAMS, depsrc=5.0, nt=600, dt=0.01, scale=1e20,
                  strike=100.0, dip=30.0, rake=70.0, zne=True, print_log=False)
    pygrt.utils.lamb(**common, rcv_points=points, output_path=output_root / "western_points")
    pygrt.utils.lamb(**common, dist=10.0, azimuth=270.0, deprcv=1.0,
                    output_path=output_root / "western_polar")
    for component in "ZNE":
        actual = read(str(output_root / "western_points" / "0000_0.00_-10.00_1.00" / f"{component}.sac"))[0].data
        expected = read(str(output_root / "western_polar" / f"{component}.sac"))[0].data
        assert np.max(np.abs(expected)) > 0
        np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=np.max(np.abs(expected))*2e-7)

    # 多点源的最早 P、S 可以来自不同源，零强度源不参与，-Ep 必须包含破裂延迟
    header = "# X-start Y-start X-fin Y-fin Kode rtlat reverse dip top bot\n1 2 3 4 5 6 7 8 9 10 11\n"
    # Kode=100 的 rake 接收断层在零滑动量时仍保留滑动角，省略剖分时取中心点
    receiver_fault = output_root / "receiver_fault.inr"
    receiver_fault.write_text(header.replace("rtlat reverse", "rake netslip") +
                              "1 8 0 10 0 100 45 0 90 .5 1.5\n" +
                              "2 8 2 10 2 100 -30 1 90 .5 1.5\n" +
                              "3 8 0 10 0 100 90 0 90 .5 1.5\n")
    for size in (None, (1, .5)):
        pygrt.utils.lamb(modelparams=MODEL_PARAMS, depsrc=5, nt=16, dt=.01, scale=1e20,
                         rcv_fault=receiver_fault, rcv_fault_size=size,
                         output_path=output_root / "receiver_faults", print_log=False)

    # 单个接收断层也使用子目录
    receiver_fault.write_text(header.replace("rtlat reverse", "rake netslip") + "1 8 0 10 0 100 45 0 90 .5 1.5\n")
    pygrt.utils.lamb(modelparams=MODEL_PARAMS, depsrc=5, nt=16, dt=.01, scale=1e20,
                     rcv_fault=receiver_fault, output_path=output_root / "receiver_fault_single", print_log=False)

    fault = output_root / "arrival_fault.inp"
    fault.write_text(header +
                     "1 4 0 6 0 400 100000 0 90 1.5 2.5 -Di+d.71\n" +
                     "2 0 0 2 0 400 100000 0 90 1.5 2.5 -Di\n" +
                     "3 8 0 10 0 400 0 0 90 1.5 2.5 -Di\n")
    receivers = output_root / "arrival_receivers.txt"
    receivers.write_text("0 10 1\n0 12 0\n")
    dt = .02
    common = dict(modelparams=MODEL_PARAMS, src_fault=fault, rcv_points=receivers,
                  dt=dt, calc_upar=True, print_log=False)
    for time_function, delays in [(None, [.72, 0]), ("i+d.13", [.14, .14])]:
        tag = "row" if time_function is None else "global"
        reference = output_root / f"arrivals_full_{tag}"
        pygrt.utils.lamb(**common, time_function=time_function, nt=700, output_path=reference)
        folders = sorted(path for path in reference.iterdir() if path.is_dir())
        for threads in [1, 4]:
            destination = output_root / f"arrivals_window_{tag}_{threads}"
            pygrt.utils.lamb(**common, time_function=time_function, nt=220, delayT0=-.2, ref_first_p=True,
                             nthreads=threads, output_path=destination)
            for ir, (east, depth) in enumerate([(10, 1), (12, 0)]):
                distances = np.hypot(east-np.array([5, 1]), 2-depth)
                arrivals = {"t0": np.min(distances / MODEL_PARAMS[0] + delays),
                            "t1": np.min(distances / MODEL_PARAMS[1] + delays)}
                begin_samples = round((arrivals["t0"]-.2)/dt)
                folder = sorted(path for path in destination.iterdir() if path.is_dir())[ir]
                for path in folder.glob("*.sac"):
                    trace = read(str(path))[0]
                    full = read(str(folders[ir] / path.name))[0]
                    np.testing.assert_allclose(trace.stats.sac.b, begin_samples*dt, rtol=2e-7)
                    for field, label in [("t0", "P"), ("t1", "S")]:
                        np.testing.assert_allclose(trace.stats.sac[field], arrivals[field], rtol=2e-7)
                        assert trace.stats.sac["k"+field] == label
                    for phase in range(2, 10):
                        assert f"t{phase}" not in trace.stats.sac
                        assert f"kt{phase}" not in trace.stats.sac
                    expected = full.data[begin_samples:begin_samples+220]
                    assert np.max(np.abs(expected)) > 0
                    np.testing.assert_allclose(trace.data, expected, rtol=3e-6, atol=np.max(np.abs(expected))*2e-7)

    # 矩形断层的初至由剖分后子源确定，移动到波形内部的时间窗仍需保留卷积历史
    rectangle = output_root / "rectangle.inp"
    rectangle.write_text(header + "1 0 0 2 0 100 .1 .2 90 1 3 -Dp/.4+d.71\n")
    common = dict(modelparams=MODEL_PARAMS, src_fault=rectangle, src_fault_size=(1, 1),
                  dist=10, azimuth=90, deprcv=0, dt=dt, calc_upar=True, print_log=False)
    reference = output_root / "rectangle_full"
    pygrt.utils.lamb(**common, nt=700, output_path=reference)
    distances = np.hypot(10-np.array([.5, 1.5, .5, 1.5]), [1.5, 1.5, 2.5, 2.5])
    first_p = np.min(distances) / MODEL_PARAMS[0] + .72
    for tag, params in [("ep", dict(delayT0=-.2, ref_first_p=True)), ("late", dict(delayT0=2.4))]:
        destination = output_root / f"rectangle_{tag}"
        pygrt.utils.lamb(**common, nt=220, nthreads=4, output_path=destination, **params)
        begin_samples = round(((first_p-.2) if tag == "ep" else 2.4)/dt)
        for path in destination.glob("*.sac"):
            if path.name == "sig.sac":
                continue
            trace = read(str(path))[0]
            full = read(str(reference / path.name))[0]
            np.testing.assert_allclose(trace.stats.sac.b, begin_samples*dt, rtol=2e-7)
            for field, speed in [("t0", MODEL_PARAMS[0]), ("t1", MODEL_PARAMS[1])]:
                np.testing.assert_allclose(trace.stats.sac[field], np.min(distances)/speed+.72, rtol=2e-7)
            expected = full.data[begin_samples:begin_samples+220]
            np.testing.assert_allclose(trace.data, expected, rtol=3e-6, atol=np.max(np.abs(expected))*2e-7)

    # 单源保留后续震相标签，每个标签都应包含显式延迟，时间函数形状不改变理论到时
    common = dict(modelparams=MODEL_PARAMS, depsrc=5, deprcv=1, dist=10, azimuth=30,
                  scale=1e20, force=(.5, -1, 2), nt=700, dt=dt, print_log=False)
    for tag, signal in [("impulse", None), ("delayed", "p/.4+d.71")]:
        pygrt.utils.lamb(**common, time_function=signal, output_path=output_root / tag)
    original = read(str(output_root / "impulse" / "Z.sac"))[0].stats.sac
    delayed = read(str(output_root / "delayed" / "Z.sac"))[0].stats.sac
    for phase in range(10):
        field = f"t{phase}"
        if field in original:
            np.testing.assert_allclose(delayed[field], original[field]+.72, rtol=2e-7)
            assert delayed["k"+field] == original["k"+field]

    # 所有子源均为零强度时，波形为零，初至及后续标签保持未定义
    silent = output_root / "silent.inp"
    silent.write_text(header + "1 0 0 2 0 400 0 0 90 1 3\n")
    pygrt.utils.lamb(modelparams=MODEL_PARAMS, src_fault=silent, rcv_points=receivers, nt=20, dt=dt,
                     delayT0=-.2, ref_first_p=True, output_path=output_root / "silent", print_log=False)
    for path in (output_root / "silent").rglob("*.sac"):
        trace = read(str(path))[0]
        assert not np.any(trace.data)
        for phase in range(10):
            assert f"t{phase}" not in trace.stats.sac
            assert f"kt{phase}" not in trace.stats.sac
