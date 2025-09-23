import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
import requests
from typing import List, Dict, Set, Tuple


_mirna_regex = re.compile(
    r'^(?P<sp>[a-z]{3})-(?P<fam>miR|mir|mirna|let)-'
    r'(?P<core>\d+[a-z]?)'
    r'(?:-(?P<extra>\d+))?'
    r'(?:-(?P<arm>[35]p))?$', re.IGNORECASE
)

server = "https://rest.ensembl.org"
endpoint = "/xrefs/id"


def check_TCGA_cols(df, output_file):
    
    tcga_columns= [col for col in df.columns if col.startswith("TCGA")]
    non_tcga_columns=[col for col in df.columns if not col.startswith("TCGA")]
    
    if "sample" in df.columns:
        filtered_df = df[["sample"] + tcga_columns]
    else:
        filtered_df = df[tcga_columns]
    
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print("file only including TCGA is saved ")
    print(f"\nNumber of non-tcga sample {len(non_tcga_columns)}")
    print("name of non-tcga columns: ")
    print(non_tcga_columns)
    

def check_TCGA_rows(df, output_file):
    
    id_column = df.columns[0]
    total_rows = df.shape[0]
    tcga_rows = df[df[id_column].str.startswith("TCGA", na=False)]
    non_tcga_count = total_rows - tcga_rows.shape[0]
    
    tcga_rows.to_csv(output_file, sep='\t', index=False)
    print("\nfile only including TCGA is saved ")
    print(f"\nNumber of non-tcga sample {non_tcga_count}")
    

def separate_normal(input_file, output_file):
    df=pd.read_csv(input_file)
    
    df["patient_id"] = df["sample"]
    df["sample_type"]=df["sample"].str.slice(13,15)
    
    normals=df[df["sample_type"]=="11"]
    primary=df[df["sample_type"] !="11"].copy()
    
    primary =primary.drop_duplicates(subset="patient_id")
    
    primary.to_csv(output_file,index=False)
    print(f"removed normal sample length: {len(normals)}")
    print(f"primary tumor sample length: {len(primary)}")
    
    
def remove_normal_samples(df):

    normals=[]
    for col in df.columns:
        if col[13:15] == "11":
            normals.append(col)
    print(f"number of removed samples is {len(normals)}")
    
    """for i in normals:  
        print(i)"""
        
    kept_col=[]
    for col in df.columns:
        if col not in normals:
            kept_col.append(col)
            
    df_cleaned= df[kept_col]
    
    return df_cleaned


def macthed_samples(df_mRNA, df_miRNA):
    mRNA_sample = df_mRNA.columns[1:]
    miRNA_sample=df_miRNA.columns[1:]
    
    mRNA_patient=mRNA_sample.str.slice(0)
    miRNA_patient=miRNA_sample.str.slice(0)
    
    common_patients =sorted(set(mRNA_patient).intersection(set(miRNA_patient)))
    print(f"number of matched patients is {len(common_patients)}")
    
    matched_mRNA_cols = []
    matched_miRNA_cols = []
    
    for patient_id in common_patients:
        for col in mRNA_sample:
            if col.startswith(patient_id):
                matched_mRNA_cols.append(col)
                break

        for col in miRNA_sample:
            if col.startswith(patient_id):
                matched_miRNA_cols.append(col)
                break
            
    matched_mRNA_df = df_mRNA[[df_mRNA.columns[0]] + matched_mRNA_cols]
    matched_miRNA_df = df_miRNA[[df_miRNA.columns[0]] + matched_miRNA_cols]

    return matched_mRNA_df, matched_miRNA_df, common_patients


def low_expression_filtering(df_miRNA: pd.DataFrame, thresold: float=0.8) ->pd.DataFrame:
    id_col = df_miRNA.columns[0]
    miRNA_ids = df_miRNA[id_col]
    expression_data = df_miRNA.iloc[:, 1:]
    
    zero_percentage = (expression_data==0).sum(axis=1)/expression_data.shape[1]
    kept_miRNA_indices = zero_percentage<thresold
    
    filtered_df_miRNA = pd.concat([miRNA_ids[kept_miRNA_indices].reset_index(drop=True),
                                   expression_data[kept_miRNA_indices].reset_index(drop=True)], axis=1)
    
    filtered_df_miRNA.columns=df_miRNA.columns
    
    removed_count = len(miRNA_ids) - len(filtered_df_miRNA)
    print(f"the number of filtering miRNAs is {removed_count}")
    
    return filtered_df_miRNA
    

def calculate_gene_wise(df, output_file):
    genes=df.iloc[:, 0]
    data=df.iloc[:,1:]
    
    stats = pd.DataFrame({
        'gene': genes,

        'min': data.min(axis=1),
        'max': data.max(axis=1),
        'q1': data.quantile(0.25, axis=1),
        'q2': data.median(axis=1),
        'q3': data.quantile(0.75, axis=1),
        'percent_zero': (data == 0).sum(axis=1) / data.shape[1] * 100
    })
    
    stats.to_excel(output_file, index=False)
    print(f"Saved gene-wise stats to: {output_file}")


def calculate_sample_wise(df, output_file):
    samples = df.columns[1:]
    data = df.iloc[:, 1:]
    
    stats = pd.DataFrame({
        "sample": samples,
        "min": data.min(),
        "max": data.max(),
        'q1': data.quantile(0.25),
        "median": data.median(),
        'q3': data.quantile(0.75),
        "zero_percentage": (data == 0).sum() / data.shape[0] * 100
    })

    stats.to_excel(output_file, index=False)
    print(f"Sample-wise statistics saved to {output_file}")


def plot_gene_expression_density(mRNA_file: str, miRNA_file:str, output_file: str=None):
    mRNA_df=pd.read_excel(mRNA_file)
    miRNA_df=pd.read_excel(miRNA_file)
    
    sns.set(style="whitegrid")
    plt.figure(figsize=(12, 6))
    
    mRNA_values = mRNA_df.values.flatten()
    miRNA_values = miRNA_df.values.flatten()  
    
    sns.kdeplot(mRNA_values, label="mRNA", fill=True,
                color="skyblue", alpha=0.5, linewidth=2)
    sns.kdeplot(miRNA_values, label="miRNA", fill=True,
                color="salmon", alpha=0.5, linewidth=2)

    plt.title("Gene-wise Expression Density (mRNA vs miRNA)", fontsize=14)
    plt.xlabel("Expression  ", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.legend()
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300)
        print(f"chart is saved as '{output_file}'")
    else:
        plt.show()

    plt.close()
    
    
def plot_full_expression_density(mRNA_file, miRNA_file, output_file=None):
    df_mRNA=pd.read_csv(mRNA_file,sep="\t")
    df_miRNA=pd.read_csv(miRNA_file,sep="\t")
        
    mRNA_val= df_mRNA.iloc[:,1:].values.flatten()
    miRNA_val=df_miRNA.iloc[:,1:].values.flatten()
        

        
    plt.figure(figsize=(10, 6))
    sns.kdeplot(mRNA_val, label="mRNA", color="skyblue", fill=True, alpha=0.4)
    sns.kdeplot(miRNA_val, label="miRNA", color="salmon", fill=True, alpha=0.4)

    plt.xlabel("Expression")
    plt.ylabel("Density")
    plt.title("Gene Expression Density (mRNA and miRNA)")
    plt.legend()
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300)
        print(f"Saved to {output_file}")
    else:
        plt.show()
            
    
def exact_match(miRNA_file, miRTarBase_file):
    df_miRNA=pd.read_csv(miRNA_file, sep="\t")
    miRNA_list=df_miRNA.iloc[:,0].unique()
    
    mirtarbase = pd.read_csv(miRTarBase_file)
    mirtarbase = mirtarbase[mirtarbase["miRNA"].str.startswith("hsa")]
    mirna_in_mirtarbase = mirtarbase["miRNA"].unique()
    
    matched = [mir for mir in miRNA_list if mir in mirna_in_mirtarbase]
    unmatched = [mir for mir in miRNA_list if mir not in mirna_in_mirtarbase]
    
    print(f"total number of miRNA: {len(miRNA_list)}")
    print(f"Eşleşen (target'i bilinen) miRNA sayisi: {len(matched)}")
    print(f"Eşleşmeyen miRNA sayisi: {len(unmatched)}")
    for elt in unmatched:
        print(elt+" ")


def _normalize_family_token(tok: str) -> str:
    # 'mir' → 'miR', 'miRNA' → 'miR', 'let' → 'let'
    t = tok.strip()
    if t.lower() in ('mir', 'mirna'):
        return 'miR'
    if t.lower() == 'let':
        return 'let'
    return t


def parse_mirna_name(name: str) -> Tuple[str, str, str, str, str]:
    """
    Dönüş: (species, family, core, extra, arm) — bulunamayanlar None döner.
    Örnekler:
      hsa-let-7a-1           -> sp=hsa, fam=let, core=7a, extra=1, arm=None
      hsa-miR-509-3          -> sp=hsa, fam=miR, core=509, extra=3, arm=None
      hsa-miR-100-5p         -> sp=hsa, fam=miR, core=100, extra=None, arm=5p
      hsa-miR-509-3-5p       -> sp=hsa, fam=miR, core=509, extra=3, arm=5p
    """
    n = name.strip().replace('—', '-').replace('–', '-')
    m = _mirna_regex.match(n)
    if not m:
        return (None, None, None, None, None)
    sp = m.group('sp').lower()
    fam = _normalize_family_token(m.group('fam'))
    core = m.group('core')
    extra = m.group('extra')
    arm = (m.group('arm') or '').lower() if m.group('arm') else None
    return (sp, fam, core, extra, arm)


def generate_mirna_candidates(name: str) -> List[str]:
    """
    - Harf suffix'lerini (a,b,...) asla atmaz.
    - Sonundaki sayılar (kopya/hairpin varyant): hem tutarak hem atarak dener.
    - Mature kolları: -3p ve -5p her taban için eklenir.
    - Çıktı: lowercase, benzersiz, deterministik sıralı.
    """
    sp, fam, core, extra, arm = parse_mirna_name(name)
    if sp is None:
        # Pars edilemeyeni olduğu gibi tek aday olarak dön (son çare)
        return [name.strip().lower()]

    # Tabanlar: (extra varsa) iki adet; yoksa bir adet
    bases = []
    fam = fam  # 'miR' veya 'let'
    if extra:
        bases.append(f"{sp}-{fam}-{core}-{extra}")
        bases.append(f"{sp}-{fam}-{core}")
    else:
        bases.append(f"{sp}-{fam}-{core}")

    # Adayları oluştur (tabanın kendisi + -3p + -5p)
    cand_set: Set[str] = set()
    for b in bases:
        cand_set.add(b)
        cand_set.add(f"{b}-3p")
        cand_set.add(f"{b}-5p")

    # Hepsini lowercase'e çevirip düzenli sırala
    cands = sorted({c.lower() for c in cand_set})
    return cands


def smart_match_mirtarbase(miRNA_file: str, miRTarBase_file: str,
                           species_filter: str = "hsa",
                           matched_csv: str = "matched_targets.csv",
                           unmatched_txt: str = "unmatched_miRNAs.txt") -> pd.DataFrame:
    """
    - miRNA dosyasındaki ilk sütunu (ID listesi) alır.
    - miRTarBase'i 'Species (miRNA) == hsa' ile filtreler.
    - Her miRNA için 3/6'lı aday isimleri üretip tam string eşleştirmesi yapar (lowercase).
    - Eşleşen tüm satırları CSV'ye yazar, eşleşmeyen ID'leri txt olarak kaydeder.
    - Dönen değer: eşleşen satırlar (DataFrame).
    """
    # miRNA listesi
    df_miRNA = pd.read_csv(miRNA_file, sep="\t")
    mirna_ids = df_miRNA.iloc[:, 0].astype(str).tolist()
    total_unique = len(set(mirna_ids))   # unique miRNAs

    # miRTarBase
    mtb = pd.read_csv(miRTarBase_file)
    if "Species (miRNA)" in mtb.columns:
        mtb = mtb[mtb["Species (miRNA)"].astype(str).str.lower() == species_filter.lower()]
    # miRTarBase miRNA sütunu zorunlu
    if "miRNA" not in mtb.columns:
        raise ValueError("miRTarBase dosyasında 'miRNA' kolonu bulunamadı.")

    # Eşleşmeyi hızlandırmak için lowercase index
    mtb = mtb.copy()
    mtb["_miRNA_lc"] = mtb["miRNA"].astype(str).str.lower()    # yardımcı sütun, hızlı ve case-insensitive

    # Küçük bir index: miRNA adı -> bu ada sahip satır indeksleri
    # aynı miRNA'den birden fazla target satırı olursa o indexler tutulut ve O(1) lookup sağlanır
    index: Dict[str, List[int]] = {}
    for i, val in enumerate(mtb["_miRNA_lc"].tolist()):
        index.setdefault(val, []).append(i)

    matched_rows: List[int] = []
    matched_by_id: Dict[str, Set[str]] = {}
    unmatched: List[str] = []

    for mid in mirna_ids:
        cands = generate_mirna_candidates(mid)
        hits_here: Set[str] = set()
        for c in cands:
            if c in index:
                hits_here.add(c)
                matched_rows.extend(index[c])
        if hits_here:
            matched_by_id[mid] = hits_here
        else:
            unmatched.append(mid)

    # unique satır indeksleri
    matched_rows = sorted(set(matched_rows))
    matched_df = mtb.iloc[matched_rows].drop(columns=["_miRNA_lc"])

    # Özet yazdır
    print(f"Toplam miRNA (unique): {total_unique}")
    print(f"Eşleşen miRNA sayısı: {len(matched_by_id)}")
    print(f"Eşleşmeyen miRNA sayısı: {len(set(unmatched))}")

    # Kaydet
    matched_df.to_csv(matched_csv, index=False)
    with open(unmatched_txt, "w", encoding="utf-8") as f:
        for u in sorted(set(unmatched)):
            f.write(u + "\n")
    print(f"Eşleşen target satırları: {matched_csv} dosyasına yazıldı.")
    print(f"Eşleşmeyen miRNA ID'leri: {unmatched_txt} dosyasına yazıldı.")

    sample_to_show = list(matched_by_id.items())[:10]
    if sample_to_show:
        print("\nÖrnek eşleşmeler (ilk 10):")
        for orig, hitset in sample_to_show:
            print(f"  {orig} -> {', '.join(sorted(hitset))}")

    return matched_df

def build_mirna_target_distribution(matched_file, output_excel):
    """
    Builds a distribution of miRNAs and their target counts from a matched targets file.

    Args:
        matched_file (str): Path to the CSV file containing matched miRNA-target pairs.
        output_excel (str): Path to the output Excel file.
    """
    try:
        matched_df = pd.read_csv(matched_file)

        if 'miRNA' not in matched_df.columns:
            raise ValueError("The 'miRNA' column is not found in the matched file.")
        
        # Count the number of unique targets for each miRNA
        target_counts = matched_df.groupby('miRNA')['Target Gene'].nunique().reset_index()
        target_counts.columns = ['miRNA', 'Target Count']
        
        # Save the result to an Excel file
        with pd.ExcelWriter(output_excel) as writer:
            target_counts.to_excel(writer, index=False, sheet_name='miRNA_Target_Counts')

        print(f"miRNA target distribution saved to '{output_excel}'")

    except FileNotFoundError:
        print(f"Error: The file '{matched_file}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")




if __name__=="__main__":

    mRNA_file=r"TCGA-BRCA.star_fpkm-uq.tsv"
    df_mRNA=pd.read_csv(mRNA_file, sep="\t")

    miRNA_file=r"TCGA-BRCA.mirna.tsv"
    df_miRNA=pd.read_csv(miRNA_file, sep="\t")

    check_TCGA_cols(df_mRNA, "cheked_mRNA")
    check_TCGA_cols(df_miRNA, "checked_miRNA")

    df_mRNA_cleaned = remove_normal_samples(df_mRNA)
    df_mRNA_cleaned.to_csv("mRNA_no_normal.tsv", sep='\t', index=False)

    df_miRNA_cleaned = remove_normal_samples(df_miRNA)
    df_miRNA_cleaned.to_csv("miRNA_no_normal.tsv", sep='\t', index=False)
    

     
    matched_mRNA_df, matched_miRNA_df, matched_patients = macthed_samples(df_mRNA_cleaned, df_miRNA_cleaned)
    
    matched_miRNA_df = low_expression_filtering(matched_miRNA_df, thresold=0.8)  
    
    matched_mRNA_df.to_csv("matched_mRNA.tsv", sep='\t', index=False)
    matched_miRNA_df.to_csv("matched_miRNA.tsv", sep='\t', index=False)

    calculate_gene_wise(matched_mRNA_df, "gene_wise_mRNA.xlsx")
    calculate_gene_wise(matched_miRNA_df, "gene_wise_miRNA.xlsx")

    calculate_sample_wise(matched_mRNA_df, "sample_wise_mRNA.xlsx")
    calculate_sample_wise(matched_miRNA_df, "sample_wise_miRNA.xlsx")

    plot_gene_expression_density("gene_wise_mRNA.xlsx", "gene_wise_miRNA.xlsx", "gene_expression_density.png")
    plot_full_expression_density("matched_mRNA.tsv", "matched_miRNA.tsv", "gene_expression_density.png")

    smart_match_mirtarbase(
    miRNA_file="matched_miRNA.tsv",         
    miRTarBase_file="MicroRNA_Target_Sites.csv", 
    species_filter="hsa",
    matched_csv="matched_targets.csv",
    unmatched_txt="unmatched_miRNAs.txt"
    )

    build_mirna_target_distribution(
        matched_file="matched_targets.csv", 
        output_excel="miRNA_target_distribution.xlsx"
    )



