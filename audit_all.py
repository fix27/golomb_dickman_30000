import sys,os,re,pickle,json,glob,math
sys.set_int_max_str_digits(0)
def cert(p):
    t=open(p,encoding="ascii").read()
    S=int(re.search(r"S=(\d+)",t).group(1)); lo=int(re.search(r"lam_lo=(\d+)",t).group(1)); hi=int(re.search(r"lam_hi=(\d+)",t).group(1))
    d=[l.strip() for l in t.splitlines() if l.strip().startswith("0.")][-1][2:]
    return S,lo,hi,d
def audit(certp, ckdir, metap, kind):
    S,lo,hi,d=cert(certp)
    if kind=="pkl":
        o=pickle.load(open(metap,"rb")); D,G,h,Ebits,KL,KR,chunk,ver=o
    else:
        o=json.load(open(metap)); pp=o["parameters"]
        D=int(pp["D"]);G=int(pp["G"]);h=float(pp["h"]);Ebits=int(pp["Ebits"]);KL=int(pp["K_L"]);KR=int(pp["K_R"]);chunk=o["chunk"]
    hn,hd=h.as_integer_ratio()
    ch=[]
    for p in glob.glob(os.path.join(ckdir,"c_*.pkl")):
        a,b,ev,od=pickle.load(open(p,"rb")); ch.append((int(a),int(b),int(ev[0]),int(ev[1]),int(od[0]),int(od[1])))
    ch.sort()
    gaps=sum(1 for i in range(len(ch)-1) if ch[i][1]!=ch[i+1][0])
    L=sum(c[2]+c[4] for c in ch); H=sum(c[3]+c[5] for c in ch)
    lo_fin=(hn*L)//hd; hi_fin=-((-hn*H)//hd)
    B=lo_fin-lo; C=hi-hi_fin
    scale=10**len(d); a=(lo*scale)>>S; b=(hi*scale)>>S
    ok=(a==b==int(d))
    print("%-46s S=%d chunks=%d nodes=%d gaps=%d B=%d C>=B=%s floorOK=%s width=B+C+%d" %
          (certp.split("\\")[-1],S,len(ch),sum(c[1]-c[0] for c in ch),gaps,B,C>=B,ok,(hi-lo)-B-C))
for args in [
 (r"30000\golomb_dickman_30000_certified.txt",r"30000\gd_rigorous_ckpt_30000",r"30000\gd_rigorous_ckpt_30000\meta.pkl","pkl"),
 (r"30000-v2\golomb_dickman_30000_certified.txt",r"30000-v2\ckpt_30000",r"30000-v2\ckpt_30000\meta_v12.json","json"),
 (r"30000-v2\golomb_dickman_2000_certified.txt",r"30000-v2\ckpt_2000",r"30000-v2\ckpt_2000\meta_v12.json","json"),
]:
    try: audit(*args)
    except Exception as e: print("FAIL",args[0],type(e).__name__,e)
