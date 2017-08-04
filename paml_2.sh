mkdir -p paml_files/two
for file_n in control/two/ctl_*; do codeml $file_n; done