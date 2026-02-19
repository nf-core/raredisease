nf-test test subworkflows/local/ --profile test,singularity --verbose --debug --update-snapshot --clean-snapshot >wflog 2>&1
nf-test test tests --profile test,singularity --verbose --debug --update-snapshot --clean-snapshot >>wflog 2>&1
