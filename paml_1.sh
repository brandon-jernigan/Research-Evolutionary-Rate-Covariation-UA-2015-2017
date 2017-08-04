mkdir -p paml_files/one
for file_n in control/one/ctl_*; do codeml $file_n; done