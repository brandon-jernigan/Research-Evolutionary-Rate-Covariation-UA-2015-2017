mkdir -p paml_files/four
for file_n in control/four/ctl_*; do codeml $file_n; done