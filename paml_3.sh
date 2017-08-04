mkdir -p paml_files/three
for file_n in control/three/ctl_*; do codeml $file_n; done