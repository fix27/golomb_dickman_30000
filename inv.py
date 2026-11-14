import os, sys, hashlib, re, glob, pickle, json
sys.set_int_max_str_digits(0)
R=os.getcwd()
def sha(p):
    return hashlib.sha256(open(p,"rb").read()).hexdigest()
def info(p):
    b=open(p,"rb").read()
    return len(b), b.count(b"\r"), sha(p)[:16]
pairs=[
 ("b084945.txt","40000/b084945.txt"),
 ("golomb_dickman_40000_certified.txt","40000/golomb_dickman_40000_certified.txt"),
 ("golomb_dickman_40000_digits.txt","40000/golomb_dickman_40000_digits.txt"),
 ("README.md","40000/README.md"),
 ("golomb_dickman_26.py","alt2-50000/golomb_dickman_26.py"),
 ("30000-v2/golomb_dickman_23.py","old_version/golomb_dickman_23.py"),
]
print("== duplicate comparison ==")
for a,b in pairs:
    pa,pb=os.path.join(R,a),os.path.join(R,b)
    if not (os.path.exists(pa) and os.path.exists(pb)):
        print("missing:",a,b); continue
    ia,ib=info(pa),info(pb)
    print("%-28s %s | %-16s %s | same=%s" % (a,ia,b.split('/')[0],ib,open(pa,'rb').read()==open(pb,'rb').read()))
print("\n== certificate headers ==")
certs=[
 "golomb_dickman_40000_certified.txt","40000/golomb_dickman_40000_certified.txt",
 "30000/golomb_dickman_30000_certified.txt","30000-v2/golomb_dickman_30000_certified.txt",
 "30000-v2/golomb_dickman_2000_certified.txt",
 "40000/smoke_60.txt","30000/smoke_60.txt","30000-v2/smoke_60.txt",
]
for c in certs:
    p=os.path.join(R,c)
    if not os.path.exists(p): print(c,"MISSING"); continue
    b=open(p,"rb").read(); t=b.decode("ascii","replace")
    m=re.search(r"S=(\d+)",t); lo=re.search(r"lam_lo=(\d+)",t); hi=re.search(r"lam_hi=(\d+)",t)
    dig=[l for l in t.splitlines() if l.startswith("0.")]
    print("%-45s bytes=%7d CR=%4d digits=%s S=%s" % (c,len(b),b.count(b"\r"),len(dig[0])-2 if dig else "?",m.group(1) if m else "?"))
print("\n== checkpoint dirs ==")
for d in ["30000/gd_rigorous_ckpt_30000","30000-v2/ckpt_30000","30000-v2/ckpt_2000","40000/gd_rigorous_ckpt_40000","alt2-50000/gd_rigorous_ckpt_50000","30000/gd_rigorous_ckpt_60","40000/gd_rigorous_ckpt_60","30000-v2/gd_rigorous_ckpt_60"]:
    p=os.path.join(R,d)
    if not os.path.isdir(p): print(d,"MISSING"); continue
    fs=os.listdir(p); cs=[f for f in fs if f.startswith("c_") and f.endswith(".pkl")]
    extra=[f for f in fs if not f.startswith("c_")]
    print("%-45s files=%4d chunks=%4d extra=%s" % (d,len(fs),len(cs),extra))
print("\n== alt2-50000 meta ==")
mp=os.path.join(R,"alt2-50000/gd_rigorous_ckpt_50000/meta_v12.json")
print(open(mp).read() if os.path.exists(mp) else "no meta")
print("\n== 30000-v2 ckpt meta? ==")
for d in ["30000-v2/ckpt_2000","30000-v2/ckpt_30000"]:
    p=os.path.join(R,d)
    print(d, [f for f in os.listdir(p) if not f.startswith("c_")][:10] if os.path.isdir(p) else "missing")
