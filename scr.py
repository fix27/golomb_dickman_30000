import re
for p in [r"40000\golomb_dickman.py", r"golomb_dickman_26.py", r"30000\golomb_dickman_rigorous_fast.py"]:
    t=open(p,encoding="utf-8",errors="replace").read()
    print("\n===",p,"===")
    for pat in [r"DIGITS\s*=\s*\d+", r"def run_certificate\([^)]*\)", r"def main\([^)]*\)", r"__main__", r"argparse", r"PROVE|prove"]:
        ms=re.findall(pat,t)
        print(pat,"->",ms[:5], "count",len(ms))
    # main block tail
    idx=t.rfind('if __name__')
    print("tail from __main__:", t[idx:idx+600].replace("\n"," | ")[:600] if idx>=0 else "none")
