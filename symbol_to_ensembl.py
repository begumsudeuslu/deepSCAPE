import pandas as pd
import requests
import time

def symbol_to_ensembl(symbols, species="hsa"):
    url="https://rest.ensembl.org/xrefs/symbol/{}/".format(species)
    headers={"Content-Type": "application/json"}

    results=[]

    for i, sym in enumerate(symbols, start=1):
        try:
            r=requests.get(url+ sym,headers=headers, timeout=10)
            r.raise_for_status()
            data=r.json()

            if data:
                ensembl_id=data[0]["id"]
            else:
                ensembl_id=None

            print(f"[{i}/{len(symbols)}] {sym} → {ensembl_id}")
            
            results.append((sym, ensembl_id))

        except Exception as e:
            print(f"Hata: {sym} için {e}")
            results.append((sym, None))
        
        time.sleep(0.1)  #küçük bekleme
    
    return pd.DataFrame(results, columns=["symbol", "ensembl_id"])

df = pd.read_csv("MicroRNA_Target_Sites.csv")
df_hsa = df[df["Species (Target Gene)"] == "hsa"].copy()

symbols = df_hsa["Target Gene"].unique().tolist()

mapping = symbol_to_ensembl(symbols, species="human")

# Orijinal tabloya ekleyelim
df = df.merge(mapping, left_on="Target Gene", right_on="symbol", how="left")
df.to_csv("mirTarBase_with_ensembl.csv", index=False)

print("Bitti! Çıktı: mirTarBase_with_ensembl.csv")