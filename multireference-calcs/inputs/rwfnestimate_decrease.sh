#!/bin/bash
rnucleus < rnucleus_input72
rwfnestimate < rwfnestimate_allTF
rmcdhf < rmcdhf_input3
rnucleus < rnucleus_input71
mv rwfn.inp rwfn.out
rwfnestimate < rwfnestimate_previous
rmcdhf < rmcdhf_input3
rnucleus < rnucleus_input70.5
mv rwfn.inp rwfn.out
rwfnestimate < rwfnestimate_previous
rmcdhf < rmcdhf_input3
# rnucleus < rnucleus_input70
# mv rwfn.inp rwfn.out
# rwfnestimate < rwfnestimate_previous
# rmcdhf < rmcdhf_input3
cp rwfn.inp rwfn.out
echo "made final rwfn.inp and rwfn.out (same file)"
