#!/usr/bin/env python3
"""fig08 + fig09 — source-spectrum validation against original data.

Re-uses the validation pipeline in benchmark/sources/:
    fig08_source_validation.pdf  all spectrum modes vs. the CosmoALEPH table
                                 and the 208-point compilation digitized from
                                 Reyna (2006) Fig. 3, with ratio panel
    fig09_charge_ratio.pdf       generator charge-ratio table vs. CosmoALEPH

The underlying CSVs are produced by benchmark/sources/extract_reyna_fig3.py
(vector-exact digitization of references/source/Reyna-B.pdf).
"""

import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SOURCES = os.path.join(REPO, "benchmark", "sources")
FIGS = os.path.join(HERE, "..", "figs")

sys.path.insert(0, SOURCES)
import plot_source_spectra as pss  # noqa: E402  (needs SOURCES on sys.path)


def main():
    comp = pss.load_compilation()
    ca = pss.load_cosmoaleph()
    pss.fig_spectrum(comp, ca)
    pss.fig_charge_ratio(ca)
    for src, dst in [("fig_source_spectrum.pdf", "fig08_source_validation.pdf"),
                     ("fig_charge_ratio.pdf", "fig09_charge_ratio.pdf")]:
        shutil.copyfile(os.path.join(SOURCES, "figures", src),
                        os.path.join(FIGS, dst))
        print(f"figs/{dst}")


if __name__ == "__main__":
    main()
