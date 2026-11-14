import os,sys,re,hashlib
sys.set_int_max_str_digits(0)
print("== normalized hash check for 40000 ==")
want={}
for ln in open(r"40000\SHA256SUMS.txt",encoding="ascii"):
    if ln.startswith("#") or not ln.strip(): continue
    h,n=ln.split(); want[n]=h
for n,h in want.items():
    p=os.path.join("40000",n)
    raw=open(p,"rb").read()
    got=hashlib.sha256(raw.replace(b"\r\n",b"\n").replace(b"\r",b"\n")).hexdigest()
    print("%-40s LF-normalized %s" % (n,"OK" if got==h else "STILL MISMATCH"))
print("\n== headers ==")
for p in [r"40000\golomb_dickman.py", r"30000\golomb_dickman_rigorous_fast.py", r"30000\golomb_dickman_digit.py"]:
    t=open(p,encoding="utf-8",errors="replace").read().splitlines()
    print("\n---",p,"lines",len(t),"---")
    for l in t[:35]: print(l[:130])
