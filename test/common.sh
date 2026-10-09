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

    cat > faults.inp <<'EOF'
  #     X-start     Y-start       X-fin       Y-fin  Kode      rt.lat     reverse   dip angle         top         bot
xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx   xxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx
  1           0           0           1           0   100         0.1         0.2          32           1           2
  2           0           0           1           0   200         0.1         0.2          47           1           2
  3           0           0           1           0   300         0.1         0.2          61           1           2
  4           0           0           1           0   400         1e6         2e6          38           1           2
  5           0           0           1           0   500         1e6         2e6          54           1           2
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
  2           0           0           1           0   100         -25        0.08          56           1           2  -Dc/0.2/0.3+d0.4
  3           0           0           1           0   100          70        0.15          48           1           2  -D0/time_function.txt+d0.3
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
    rm -f faults_rupture.inr faults_partial.inr faults_negative_delay.inr
}

# Check the exit status only; do not match diagnostic messages
expect_fail() {
    if "$@"; then
        echo "ERROR: invalid arguments were accepted: $*" >&2
        exit 1
    fi
}
