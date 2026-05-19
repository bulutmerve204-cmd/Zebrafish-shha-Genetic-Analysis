from Bio import SeqIO
import pandas as pd

# 1. We load the data
file =' ncbi_dataset '

fasta_yolu = "rna.fna"
#2.Veriyi Oku ve Listeye Aktar
gen_verileri = []

for kayit in SeqIO.parse(fasta_yolu, "fasta"):
    dizi = str(kayit.seq)
    uzunluk= len(dizi)

    g_sayisi = dizi.count("G")
    c_sayisi =dizi.count("C")
    gc_orani = (g_sayisi + c_sayisi) / uzunluk * 100
    gen_verileri.append({
        "Erisim_No" : kayit.id,
        "Tanim":kayit.description [:50] +"...",
        "Uzunluk_bp": uzunluk,
        "GC_Orani_%": round(gc_orani,2)
})
    df = pd.DataFrame(gen_verileri)
print(df)
import matplotlib.pyplot as plt

df.plot(kind="bar" , x = "Erisim_No" , y="GC_Orani_%", color="skyblue", title="GC_Oranlari")
plt.show()
