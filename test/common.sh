# Create shared example inputs in the current module directory
# ===================================================================================
create_test_files() {
    cat > dists <<'EOF'
 2.0
 5.0
10.0
EOF

    cat > depsrc <<'EOF'
1.0
2.0
EOF

    cat > deprcv <<'EOF'
0.0
1.0
EOF

    # The first numeric column in every sCoulomb fault row is a placeholder, always 1
    cat > faults.inp <<'EOF'
  #     X-start     Y-start       X-fin       Y-fin  Kode      rt.lat     reverse   dip angle         top         bot
xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx   xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx
  1           0           0           1           0   100         0.1         0.2          32           1           2
  1           0           0           1           0   200         0.1         0.2          47           1           2
  1           0           0           1           0   300         0.1         0.2          61           1           2
  1           0           0           1           0   400         1e6         2e6          38           1           2
  1           0           0           1           0   500         1e6         2e6          54           1           2
EOF

    awk 'NR <= 2 || $6 == 100 || $6 == 200 || $6 == 300' faults.inp > faults_rectangles.inp
    awk 'NR <= 3' faults.inp > fault_single.inp

    # A surface-reaching thrust fault for subdivision and synthesis
    cat > faults_surface.inr <<'EOF'
  #     X-start     Y-start       X-fin       Y-fin  Kode        rake     netslip   dip angle         top                bot
xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx   xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxxxxxxxxx
  1         -10           0          10           0   100          90         0.1          45           0  7.071067811865475
EOF

    # A shuffled plane with three strike columns and two dip rows, rounded to six decimal places
    cat > fault_plane.inr <<'EOF'
  #     X-start     Y-start       X-fin       Y-fin  Kode        rake     netslip   dip angle         top         bot
xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx   xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx
  1          11          20          12          20   100          15         0.2          60           2    2.866025  -Dc/1/1+d0.2
  1          10        19.5          11        19.5   100         -25         0.3          60    2.866025    3.732051  -Di+d0.3
  1          10          20          11          20   100          45           0          60           2    2.866025
  1          11        19.5          12        19.5   100          70         0.4          60    2.866025    3.732051  -D0/unused.txt+d0.4
  1          12          20          13          20   100          30         0.5          60           2    2.866025  -Di+d0.5
  1          12        19.5          13        19.5   100         -60         0.6          60    2.866025    3.732051  -Di+d0.6
EOF

    cat > faults.inr <<'EOF'
  #     X-start     Y-start       X-fin       Y-fin  Kode        rake     netslip   dip angle         top         bot
xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx   xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx
  1           0           0           1           0   100          45         0.1          43           1           2
EOF

    cat > rcv_faults.inr <<'EOF'
  #     X-start     Y-start       X-fin       Y-fin  Kode        rake     netslip   dip angle         top         bot
xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx   xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx
  1           8           0           9           0   100          45           0          57         0.2         0.8
EOF

    cat > faults_rupture.inr <<'EOF'
  #     X-start     Y-start       X-fin       Y-fin  Kode        rake     netslip   dip angle         top         bot  time_function
xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx   xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxxxxx
  1           0           0           1           0   100          35        0.12          39           1           2  -Dp/0.5+d0.2
  1           0           0           1           0   100         -25        0.08          56           1           2  -Dc/0.2/0.3+d0.4
  1           0           0           1           0   100          70        0.15          48           1           2  -D0/time_function.txt+d0.3
EOF

    sed '4s/  *-D.*//' faults_rupture.inr > faults_partial.inr
    sed 's/+d0.2/+d-0.2/' faults_rupture.inr > faults_negative_delay.inr

    cat > rcv_points.txt <<'EOF'
# north  east  depth
      0     5      0
      0    10      1
EOF

    cat > rcv_geometry.txt <<'EOF'
# north  east  depth  strike  dip  rake
      0     5      0      10   20    30
      0    10      1      40   50    60
EOF
}
# ===================================================================================

# Remove the example inputs created above
remove_test_files() {
    rm -f dists depsrc deprcv faults.inp faults.inr rcv_faults.inr rcv_points.txt rcv_geometry.txt
    rm -f faults_rupture.inr faults_partial.inr faults_negative_delay.inr faults_rectangles.inp
    rm -f fault_single.inp fault_plane.inr faults_surface.inr
}

# Check the exit status only; do not match diagnostic messages
expect_fail() {
    if "$@"; then
        echo "ERROR: invalid arguments were accepted: $*" >&2
        exit 1
    fi
}
