import shutil
import subprocess
from pathlib import Path

import pygrt

dist = 10.0
depsrc = 2.0
deprcv = 3.0
nt = 600
dt = 0.02
modname = "../milrow"
az = 22.0

pymod = pygrt.PyModel1D(grn="GRN", modelpath=modname)
pymod.greenfn(depsrc=depsrc, deprcv=deprcv, dists=dist, nt=nt, dt=dt, calc_upar=True)

pymod.syn(dist=dist, azimuth=az, scale=1e20, output_path="syn", calc_upar=True)
pygrt.utils.strain("syn")
pygrt.utils.stress("syn")
pygrt.utils.rotation("syn")

pymod.syn(dist=dist, azimuth=az, scale=1e20, output_path="syn_zne", zne=True, calc_upar=True)
pygrt.utils.strain("syn_zne")
pygrt.utils.stress("syn_zne")
pygrt.utils.rotation("syn_zne")

# 检查 -G 的两种写法和旧接口警告，张量结果须逐字节一致
for directory in ["syn", "syn_zne"]:
    for module in ["strain", "stress", "rotation"]:
        saved = {path: path.read_bytes() for path in Path(directory).glob(f"{module}_*.sac")}
        for args, legacy in [([f"-G{directory}"], False), (["-G", directory], False), ([directory], True)]:
            result = subprocess.run(["grt", module, *args], capture_output=True, text=True)
            assert result.returncode == 0, result.stdout+result.stderr
            if legacy:
                assert "[WARNING]" in result.stdout and "deprecated" in result.stdout and "-G<syn_dir>" in result.stdout
            else:
                assert "[WARNING]" not in result.stdout
            for path, content in saved.items():
                assert path.read_bytes() == content

        # 缺少目录、未知选项、混用两种目录参数和附加参数应在改写结果前失败
        for args in [[], ["-G"], ["-x"], [directory, f"-G{directory}"], [f"-G{directory}", "unused_directory"]]:
            result = subprocess.run(["grt", module, *args], capture_output=True, text=True)
            assert result.returncode != 0, args
            for path, content in saved.items():
                assert path.read_bytes() == content

# -------------------- 静态应变 / 应力 / 旋转 --------------------
pymod_s = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=modname)
pymod_s.static_greenfn(depsrc=2.0, deprcv=0.0, norths=[-3.0, 3.0, 1.0], easts=[-2.0, 2.0, 1.0], calc_upar=True)
pymod_s.static_syn(scale=1e20, output_path="stsyn.nc", calc_upar=True)
pygrt.utils.static_strain("stsyn.nc")
pygrt.utils.static_stress("stsyn.nc")
pygrt.utils.static_rotation("stsyn.nc")

pymod_s.static_syn(scale=1e20, output_path="stsyn_zne.nc", zne=True, calc_upar=True)
pygrt.utils.static_strain("stsyn_zne.nc")
pygrt.utils.static_stress("stsyn_zne.nc")
pygrt.utils.static_rotation("stsyn_zne.nc")

rcv = Path("rcv_pts.txt")
rcv.write_text("# north east depth (km)\n0 0 0\n1 2 0\n-1 1 0\n")
pymod_s.static_syn(scale=1e20, output_path="stsyn_q.nc", rcv_points=rcv, calc_upar=True)
pygrt.utils.static_strain("stsyn_q.nc")
pygrt.utils.static_stress("stsyn_q.nc")
pygrt.utils.static_rotation("stsyn_q.nc")

for name in ["GRN", "syn", "syn_zne", "stgrn.nc", "stsyn.nc", "stsyn_zne.nc", "stsyn_q.nc", "rcv_pts.txt"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
