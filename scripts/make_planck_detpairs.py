#!/usr/bin/env python3
import numpy as np

def list_planck(detset, good=True, subset=None, extend_857=True, extend_545=False):
    detectors = []
    if subset is None:
        subset = 0

    if detset in (30, "30", "030", "30GHz", "030GHz", "30A", "030A", "30B", "030B"):
        horns = range(27, 29)
        instrument = "LFI"
    elif detset in (44, "44", "044", "44GHz", "044GHz", "44A", "044A", "44B", "044B"):
        horns = range(24, 27)
        instrument = "LFI"
    elif detset in (70, "70", "070", "70GHz", "070GHz"):
        horns = range(18, 24)
        if subset == 1:
            horns = [18, 23]
        elif subset == 2:
            horns = [19, 22]
        elif subset == 3:
            horns = [20, 21]
        instrument = "LFI"
    elif detset in ["70A", "070A"]:
        horns = [18, 20, 23]
        instrument = "LFI"
    elif detset in ["70B", "070B"]:
        horns = [19, 21, 22]
        instrument = "LFI"
    elif isinstance(detset, str) and detset.upper() == "LFI":
        detectors.extend(list_planck(30, good=good))
        detectors.extend(list_planck(44, good=good))
        detectors.extend(list_planck(70, good=good))
        return detectors
    elif detset in (100, "100", "100GHz"):
        psb_horns = range(1, 5)
        swb_horns = []
        instrument = "HFI"
        freq = "100-"
    elif detset == "100A":
        psb_horns = [1, 4]
        swb_horns = []
        instrument = "HFI"
        freq = "100-"
    elif detset == "100B":
        psb_horns = [2, 3]
        swb_horns = []
        instrument = "HFI"
        freq = "100-"
    elif detset in (143, "143", "143GHz"):
        psb_horns = np.arange(1, 5)
        if good:
            swb_horns = range(5, 8)
        else:
            swb_horns = range(5, 9)
        if subset == 1:
            psb_horns, swb_horns = [1, 3], []
        elif subset == 2:
            psb_horns, swb_horns = [2, 4], []
        elif subset == 3:
            psb_horns, swb_horns = [], [5, 6, 7]
        instrument = "HFI"
        freq = "143-"
    elif detset == "143A":
        psb_horns = [1, 3]
        swb_horns = [5, 7]
        instrument = "HFI"
        freq = "143-"
    elif detset == "143B":
        psb_horns = [2, 4]
        swb_horns = [6]
        instrument = "HFI"
        freq = "143-"
    elif detset in (217, "217", "217GHz"):
        psb_horns = np.arange(5, 9)
        swb_horns = np.arange(1, 5)
        if subset == 1:
            psb_horns, swb_horns = [5, 7], []
        elif subset == 2:
            psb_horns, swb_horns = [6, 8], []
        elif subset == 3:
            psb_horns, swb_horns = [], [1, 2, 3, 4]
        instrument = "HFI"
        freq = "217-"
    elif detset == "217A":
        psb_horns = [5, 7]
        swb_horns = [1, 3]
        instrument = "HFI"
        freq = "217-"
    elif detset == "217B":
        psb_horns = [6, 8]
        swb_horns = [2, 4]
        instrument = "HFI"
        freq = "217-"
    elif detset in (353, "353", "353GHz"):
        psb_horns = np.arange(3, 7)
        swb_horns = [1, 2, 7, 8]
        if subset == 1:
            psb_horns, swb_horns = [3, 5], []
        elif subset == 2:
            psb_horns, swb_horns = [4, 6], []
        elif subset == 3:
            psb_horns, swb_horns = [], [1, 2, 7, 8]
        instrument = "HFI"
        freq = "353-"
    elif detset == "353A":
        psb_horns = [3, 5]
        swb_horns = [1, 7]
        instrument = "HFI"
        freq = "353-"
    elif detset == "353B":
        psb_horns = [4, 6]
        swb_horns = [2, 8]
        instrument = "HFI"
        freq = "353-"
    elif detset in (545, "545", "545GHz"):
        psb_horns = []
        if good and not extend_545:
            swb_horns = [1, 2, 4]
        else:
            swb_horns = np.arange(1, 5)
        instrument = "HFI"
        freq = "545-"
    elif detset == "545A":
        psb_horns = []
        swb_horns = [1]
        instrument = "HFI"
        freq = "545-"
    elif detset == "545B":
        psb_horns = []
        swb_horns = [2, 4]
        instrument = "HFI"
        freq = "545-"
    elif detset in (857, "857", "857GHz"):
        psb_horns = []
        if good and not extend_857:
            swb_horns = [1, 2, 3]
        else:
            swb_horns = np.arange(1, 5)
        instrument = "HFI"
        freq = "857-"
    elif detset == "857A":
        psb_horns = []
        swb_horns = [1, 3]
        instrument = "HFI"
        freq = "857-"
    elif detset == "857B":
        psb_horns = []
        swb_horns = [2, 4]
        instrument = "HFI"
        freq = "857-"
    elif isinstance(detset, str) and detset.upper() == "HFI":
        detectors.extend(list_planck(100, good=good, extend_857=extend_857))
        detectors.extend(list_planck(143, good=good, extend_857=extend_857))
        detectors.extend(list_planck(217, good=good, extend_857=extend_857))
        detectors.extend(list_planck(353, good=good, extend_857=extend_857))
        detectors.extend(list_planck(545, good=good, extend_857=extend_857))
        detectors.extend(list_planck(857, good=good, extend_857=extend_857))
        return detectors
    elif isinstance(detset, str) and detset.upper() == "PLANCK":
        detectors.extend(list_planck("LFI", good=good, extend_857=extend_857))
        detectors.extend(list_planck("HFI", good=good, extend_857=extend_857))
        return detectors
    else:
        # single detectors and horns
        lfidets = list_planck("LFI")
        hfidets = list_planck("HFI")
        if detset in lfidets or detset in hfidets:
            return [detset]
        if detset + "M" in lfidets:
            return [detset + "M", detset + "S"]
        if detset + "a" in hfidets:
            return [detset + "a", detset + "b"]
        return -1

    if instrument == "LFI":
        for horn in horns:
            for arm in ["S", "M"]:
                detectors.append("LFI" + str(horn) + arm)
    elif instrument == "HFI":
        for horn in psb_horns:
            for arm in ["a", "b"]:
                detectors.append(freq + str(horn) + arm)
        for horn in swb_horns:
            detectors.append(freq + str(horn))

    return detectors


def build_detpairs():
    freqs = [100, 143, 217, 353]
    detsets = [f"{freq:03}{suffix}"
               for suffix in ["GHz", "A", "B"]
               for freq in freqs]

    detsetpairs = []

    # Full frequency and detector-set auto and cross spectra
    for detset1 in detsets:
        for detset2 in detsets:
            # No cross between full-frequency and A/B subsets
            if detset1.endswith("GHz") and detset2[-1] in "AB":
                continue
            if detset2.endswith("GHz") and detset1[-1] in "AB":
                continue
            detsetpairs.append((detset1, detset2))

    # Single detector and single horn auto spectra
    for det in list_planck("Planck"):
        detsetpairs.append((det, det))
        if det.endswith("a") or det.endswith("M"):
            horn = det[:-1]
            detsetpairs.append((horn, horn))

    return detsetpairs


if __name__ == "__main__":
    pairs = build_detpairs()
    print("detpairs:")
    for d1, d2 in pairs:
        print(f"  - ['{d1}', '{d2}']")
