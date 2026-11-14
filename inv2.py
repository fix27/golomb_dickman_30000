import os,sys,re,hashlib,difflib
sys.set_int_max_str_digits(0)
R=os.getcwd()
def digits(p):
    t=open(p,encoding="ascii").read()
    lines=[l.strip() for l in t.splitlines() if l.strip()]
    for l in reversed(lines):
        if l.startswith("0."): return l[2:]
    return "".join(lines)
D40=digits(r"40000\golomb_dickman_40000_digits.txt")
print("len 40000 digits:",len(D40))
for p in ["30000\\golomb_dickman_30000_digits.txt","30000-v2\\golomb_dickman_30000_certified.txt","30000-v2\\golomb_dickman_2000_certified.txt","30000\\golomb_dickman_30000_certified.txt"]:
    d=digits(p); n=min(len(d),len(D40))
    pref=next((i for i in range(n) if d[i]!=D40[i]),n)
    print("%-58s len=%5d prefix-match=%5d" % (p,len(d),pref))
for p in ["40000\\smoke_60.txt","30000\\smoke_60.txt","30000-v2\\smoke_60.txt"]:
    d=digits(p); print("%-40s ==prefix60: %s  %s" % (p,d==D40[:60],d[:30]))
print("\n== sha256 of 40000 bundle vs SHA256SUMS ==")
want={}
for ln in open(r"40000\SHA256SUMS.txt",encoding="ascii"):
    if ln.startswith("#") or not ln.strip(): continue
    h,n=ln.split(); want[n]=h
for n,h in want.items():
    p=os.path.join(R,"40000",n)
    got=hashlib.sha256(open(p,"rb").read()).hexdigest()
    print("%-40s %s" % (n,"OK" if got==h else "MISMATCH"))
print("\n== script normalized comparison ==")
def norm(p): return open(p,"rb").read().replace(b"\r\n",b"\n")
a=norm(r"golomb_dickman_26.py"); b=norm(r"alt2-50000\golomb_dickman_26.py")
print("root26 vs alt2 norm equal:",a==b,"sizes",len(a),len(b))
if a!=b:
    da=a.decode("utf-8","replace").splitlines(); db=b.decode("utf-8","replace").splitlines()
    sm=difflib.SequenceMatcher(None,da,db); print("ratio",round(sm.ratio(),4))
    for tag,i1,i2,j1,j2 in sm.get_opcodes():
        if tag!="equal":
            print(tag,(i1+1,i2),(j1+1,j2))
            for l in da[i1:i2][:3]: print("  A:",l[:110])
            for l in db[j1:j2][:3]: print("  B:",l[:110])
