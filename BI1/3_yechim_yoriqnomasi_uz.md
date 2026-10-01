# pUC19 loyihasi: yechim yo'riqnomasi (buyruqlar, skript, P2P javoblari)

> **Muhim eslatma.** Bu yo'riqnoma har bir topshiriq uchun to'liq ish tartibini, tayyor buyruqlarni, ishlaydigan skriptni va P2P savollariga javob asoslarini beradi. Lekin 1-va 2-topshiriqlardagi **aniq FastQC natijalari** (qaysi modul qizil/to'q sariq, qaysi adapter topildi) yuklab olingan **real ma'lumotga** bog'liq — ularni men ko'rmadim va uydirmayman. Bu natijalarni o'zing FastQC hisobotidan olib, pastdagi jadvallarga yozasan. Ko'chirilgan yechim P2P'da ishlamaydi: hamma narsani o'zing tushuntira olishing kerak.

## Mundarija

1. [Tayyorgarlik](#1-tayyorgarlik)
2. [1-topshiriq](#2-1-topshiriq-fastqc-bilan-sifatni-baholash)
3. [2-topshiriq](#3-2-topshiriq-trimmomatic-bilan-tozalash)
4. [3-topshiriq](#4-3-topshiriq-gc-tarkib)
5. [Topshirish (Git)](#5-topshirish-git)
6. [Yakuniy tekshiruv ro'yxati](#6-yakuniy-tekshiruv-royxati)

---

## 1. Tayyorgarlik

### 1.1. Muhitni tekshirish

```bash
python3 --version        # 3.8+
java -version            # Trimmomatic uchun
fastqc --version
pip install biopython
python3 -c "import Bio; print(Bio.__version__)"
```

Hammasini conda bilan o'rnatish (qulay yo'l):

```bash
conda create -n puc19 -c bioconda -c conda-forge python=3.11 fastqc trimmomatic sra-tools biopython
conda activate puc19
```

### 1.2. Repozitoriy va papkalar

README talabiga ko'ra:

```bash
git clone <GitLab-repo-URL>
cd <repo>
git checkout -b develop
mkdir -p src/project1/{results,data,script}
```

Natijada struktura:

```
src/
└── project1/
    ├── data/       ← .fastq, adapter fayllari, fasta1, fasta2
    ├── results/    ← FastQC / Trimmomatic natijalari (1.1.txt, 2.1.txt ...)
    └── script/     ← gc_content.py
```

---

## 2. 1-topshiriq: FastQC bilan sifatni baholash

### 2.1. Ma'lumotni yuklash (1.2-band)

1. https://www.ncbi.nlm.nih.gov/sra/SRX1684458[accn] sahifasini och.
2. Sahifadan quyidagilarni **yozib ol** (keyin kerak bo'ladi): **Platform / Instrument**, **Layout** (SINGLE yoki PAIRED), **read uzunligi**. Bu Trimmomatic rejimini (`SE`/`PE`) va adapter faylini tanlashga ta'sir qiladi.
3. **Run** identifikatorini (`SRR...`) top.
4. Yuklash:

```bash
cd src/project1/data
prefetch SRRxxxxxxx
fasterq-dump SRRxxxxxxx -O .
ls -lh        # .fastq fayl(lar) paydo bo'lganini tekshir
```

`SRRxxxxxxx` — sahifada ko'rgan haqiqiy Run ID. Agar **PAIRED** bo'lsa, `_1.fastq` va `_2.fastq` hosil bo'ladi; loyiha tavsifi bitta fayl haqida gapirgani uchun, ehtimol single-end — lekin sahifadan o'zing tasdiqla.

Tez tekshiruv:

```bash
head -n 4 SRRxxxxxxx.fastq                       # 4 qator = 1 read
echo $(( $(wc -l < SRRxxxxxxx.fastq) / 4 ))      # read'lar soni
```

### 2.2. FastQC (1.3-band)

```bash
cd ..   # project1/ ga qaytish
fastqc data/SRRxxxxxxx.fastq -o results/
```

Fayl nomini topshiriq talabi bo'yicha qo'y (misol: `1.1.txt`). Sodda yo'l — `summary.txt`ni nomlab nusxalash:

```bash
cd results
unzip SRRxxxxxxx_fastqc.zip
cp SRRxxxxxxx_fastqc/summary.txt 1.1.txt
cat 1.1.txt
```

Natija ko'rinishi (har qatorda status, modul, fayl):

```
PASS    Basic Statistics                SRRxxxxxxx.fastq
WARN    Per base sequence quality       SRRxxxxxxx.fastq
...
```

### 2.3. Qizil va to'q sariq statuslar (1.4-band)

Terminaldan oling:

```bash
grep "^FAIL" 1.1.txt     # qizil (kritik)
grep "^WARN" 1.1.txt     # to'q sariq (e'tibor talab qiladi)
```

Va quyidagi jadvalga o'zing yoz (men bo'sh qoldirdim — to'ldirish sening vazifang):

| Modul | Status (PASS/WARN/FAIL) |
|-------|-------------------------|
| Basic Statistics | |
| Per base sequence quality | |
| Per tile sequence quality | |
| Per sequence quality scores | |
| Per base sequence content | |
| Per sequence GC content | |
| Per base N content | |
| Sequence length distribution | |
| Sequence duplication levels | |
| Overrepresented sequences | |
| Adapter content | |

Natijani `results/1.4.txt` ga saqla.

### 2.4. Adapterlar (1.5-band)

1. HTML hisobotda **Adapter Content** grafigini och: legendada qaysi adapter turlari ko'rsatilgan (masalan, *Illumina Universal Adapter*, *Illumina Small RNA Adapter*, *Nextera Transposase Sequence*).
2. **Overrepresented sequences** jadvalida **Possible Source** ustunini o'qi.
3. Nomlarni `results/1.5.txt` ga yoz.

Qo'shimcha: jadvaldagi eng ko'p uchraydigan ketma-ketlikni terminalda adapter fayllari bilan solishtirish mumkin:

```bash
grep -i "<overrepresented ketma-ketlikning boshi>" adapters/*.fa
```

(Trimmomatic'ning `adapters/` papkasidagi `.fa` fayllarida shu adapterlar ketma-ketligi bor.)

---

## 3. 2-topshiriq: Trimmomatic bilan tozalash

### 3.1. Adapter faylini tanlash (2.1-band)

Trimmomatic o'rnatilgan joyda `adapters/` papkasi bor (conda'da: `$CONDA_PREFIX/share/trimmomatic*/adapters/`). Topish:

```bash
find $CONDA_PREFIX -name "TruSeq3-SE.fa" 2>/dev/null
```

Mos faylni `data/` ga nusxala:

```bash
cp <topilgan-yo'l>/TruSeq3-SE.fa data/
```

**Tanlash mantig'i:** 1.5-bandda topilgan adapter nomi + NCBI sahifasidagi instrument/layout. Misol: Illumina HiSeq/MiSeq + single-end → `TruSeq3-SE.fa`; Nextera adapteri topilsa → `NexteraPE-PE.fa`; eski GAII → `TruSeq2-SE.fa`.

### 3.2. Birinchi urinish (2.2-band)

```bash
trimmomatic SE -phred33 \
    data/SRRxxxxxxx.fastq results/2.2_trimmed.fastq \
    ILLUMINACLIP:data/TruSeq3-SE.fa:2:30:10
```

(Agar `trimmomatic` buyrug'i yo'q bo'lsa: `java -jar /yo'l/trimmomatic-0.39.jar SE ...`)

Konsol chiqishida ko'rinadi: `Input Reads`, `Surviving`, `Dropped`. Buni `results/2.2.txt` ga saqla:

```bash
trimmomatic SE -phred33 data/SRRxxxxxxx.fastq results/2.2_trimmed.fastq \
    ILLUMINACLIP:data/TruSeq3-SE.fa:2:30:10 2> results/2.2.txt
```

Bu bosqichda faqat ILLUMINACLIP ishlatilgani uchun adapterlarning yo'qolishi toza ko'rinadi (sifat bo'yicha kesish yo'q).

### 3.3. FastQC qayta (2.3-band)

```bash
fastqc results/2.2_trimmed.fastq -o results/
cd results
unzip -o 2.2_trimmed_fastqc.zip
cp 2.2_trimmed_fastqc/summary.txt 2.3.txt
```

### 3.4. Adapterlar tozalandimi? (2.4-band)

`2.3.txt` ichida **Adapter Content** qatoriga qara:

- `PASS` — tozalandi. Fayl to'g'ri tanlangan.
- `WARN`/`FAIL` — tozalanmadi. HTML'da grafik hamon 0% dan baland bo'lsa, boshqa adapter faylini sinab ko'r va qadamlarni takrorla.

Eng yaxshi natija bergan fayl `data/` da qolsin, boshqalarini o'chir (README: ortiqcha fayl bo'lmasin).

> Agar bir nechta faylni sinasang, har bir urinishni alohida nom bilan saqla (`2.2_try1`, `2.2_try2`), va P2P'da nima uchun aynan oxirgisini tanlaganingni tushuntir.

### 3.5. Qolgan kritik muammolar (2.5-band)

1.3-bandning jadvalini shu yangi hisobot (`2.3.txt`) uchun qayta to'ldir va `results/2.5.txt` ga saqla. Qaysi modullar yaxshilandi, qaysilari o'zgarmadi, qaysilari yomonlashdi — solishtir.

### 3.6. Uzunlik taqsimoti (2.6-band, P2P savol)

**Savol:** «Nima uchun read'lar uzunligi taqsimoti o'zgardi?»

**Javob asosi:**

> Xom ma'lumotda barcha read'lar bir xil uzunlikda (sekvenator har bir read'ni qat'iy sikllar soni davomida o'qiydi), shuning uchun *Sequence length distribution* grafigida bitta tik cho'qqi bor. Trimmomatic har bir read'ni **alohida** qayta ishlaydi: ILLUMINACLIP adapter topilgan read'larning adapter qismini kesib tashlaydi (adapter qayerdan boshlansa, o'sha joydan oxirigacha), SLIDINGWINDOW/TRAILING esa past sifatli oxirlarni kesadi. Adapteri yo'q va sifati yuqori read'lar o'zgarmaydi, boshqalari esa turli darajada qisqaradi. Natijada read'lar turli uzunlikda bo'lib qoladi va grafik bitta cho'qqidan keng taqsimotga aylanadi. Bu texnik xato emas, balki tozalashning kutilgan oqibati: ma'lumot sifati uzunlik hisobiga yaxshilandi. FastQC bu modulga WARN berishi mumkin, chunki u bir xil uzunlikni «normal» deb hisoblaydi.

**O'zing tekshir:** `2.3.txt`/HTML'da haqiqatan ham uzunlik taqsimoti kengayganini ko'r. Agar kengaymagan bo'lsa (masalan, adapter deyarli topilmagan), javobingni shunga moslash.

### 3.7. Barcha funksiyalar bilan standart parametrlar (2.7-band)

Manualdagi mashhur standart ko'rsatmalar to'plami (Trimmomatic hujjatidagi namuna asosida):

```bash
trimmomatic SE -phred33 \
    data/SRRxxxxxxx.fastq results/2.7_trimmed.fastq \
    ILLUMINACLIP:data/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36 2> results/2.7.txt
```

Keyin FastQC:

```bash
fastqc results/2.7_trimmed.fastq -o results/
cd results
unzip -o 2.7_trimmed_fastqc.zip
cp 2.7_trimmed_fastqc/summary.txt 2.7_fastqc.txt
```

**Parametrlar mazmuni:**

| Parametr | Qiymat | Sababi |
|----------|--------|--------|
| ILLUMINACLIP | `2:30:10` | Standart sozlama; 2 ta mismatch ruxsat, adapterni aniq topish ostonalari |
| LEADING | 3 | Boshidagi juda past sifatli (Q<3) asoslarni kesadi |
| TRAILING | 3 | Oxiridagi juda past sifatli asoslarni kesadi |
| SLIDINGWINDOW | `4:15` | 4 nt oynada o'rtacha Q<15 bo'lsa kesadi (xato ehtimoli ~3%) |
| MINLEN | 36 | 36 nt dan qisqa read'lar alignment uchun yaroqsiz, tashlanadi |

> `MINLEN:36` qiymati read uzunligingga bog'liq: agar read'lar 36 nt dan qisqa bo'lsa, hammasi tashlanadi! Avval read uzunligini 1.2-bandda yozib olganing qiymatdan tekshir va MINLEN'ni uning ~1/2–2/3 qismiga qo'y.

### 3.8. Parametrlar haqida P2P javobi (2.8-band)

**Savol:** «Tanlagan parametrlaringni tasvirla. Nima uchun aynan shunday tanlov? Hisobot yaxshilandimi (ha/yo'q) va nima uchun?»

**Javob tuzilmasi (har bir bandni o'z natijalaring bilan to'ldir):**

1. **Tanlangan parametrlar:** 3.7-dagi jadval + qaysi adapter fayli va nima uchun.
2. **Tanlov asoslari:** FastQC dastlabki hisobotida nima ko'rdim (masalan, oxirida sifat pasayishi → TRAILING/SLIDINGWINDOW; adapter → ILLUMINACLIP; qisqa «qoldiqlar» → MINLEN).
3. **Natija (ha/yo'q):** `1.1.txt` va `2.7_fastqc.txt`ni modul bo'yicha solishtir. Masalan, *Per base sequence quality* qanday o'zgardi, *Adapter content*, *Sequence length distribution*.
4. **Nima uchun shunday:** yaxshilangan modullar uchun sabab (kesilgan sifatsiz uchlar). Yaxshilanmagan modul uchun ham sabab bo'lishi kerak: masalan, *Per base sequence content* tozalash bilan tuzalmaydi, chunki bu kutubxona tayyorlash artefakti (boshidagi bir necha pozitsiyalarda tasodifiy bo'lmagan nukleotid ulushi); *Sequence duplication levels* esa pUC19 kabi kichik genomda juda yuqori coverage tufayli yuqori bo'lishi tabiiy.
5. **Savdo-sotiq (trade-off):** qat'iy parametrlar sifatni oshiradi, lekin read'lar sonini kamaytiradi. `Surviving` foizini (`2.7.txt`) keltir.

---

## 4. 3-topshiriq: GC-tarkib

### 4.1. Fayllarni yuklash

NCBI sahifalarida: **Send to → Complete Record → File → Format: FASTA → Create File**. Yoki terminalda (Entrez orqali):

```bash
cd src/project1/data
```

Eng ishonchli yo'l — brauzerdan yuklab, quyidagicha nomlash:

```
data/fasta1.fasta    ← https://www.ncbi.nlm.nih.gov/nuccore/6691170
data/fasta2.fasta    ← https://www.ncbi.nlm.nih.gov/nuccore/MT856194.1
```

Tekshiruv:

```bash
head -n 2 data/fasta1.fasta
grep -c "^>" data/fasta1.fasta data/fasta2.fasta     # har birida yozuvlar soni
```

### 4.2. Skript (3.3-band)

Quyidagi skriptni `script/gc_content.py` ga saqla. Skript sinab ko'rilgan (kichik misolda: `GGCCAATT` + `acgtN` → natijalar to'g'ri chiqdi).

```python
#!/usr/bin/env python3
"""FASTA fayl(lar)idagi GC-tarkibni (%) hisoblaydi.

Ishlatilishi:
    python3 gc_content.py ../data/fasta1.fasta ../data/fasta2.fasta
"""
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Biopython o'rnatilmagan: pip install biopython")


def gc_stats(sequence: str) -> dict:
    """Ketma-ketlik uchun GC statistikasini qaytaradi."""
    seq = sequence.upper()
    g, c = seq.count("G"), seq.count("C")
    a, t = seq.count("A"), seq.count("T")
    acgt = a + t + g + c
    total = len(seq)
    return {
        "length": total,
        "A": a, "T": t, "G": g, "C": c,
        "other": total - acgt,  # N va boshqa IUPAC belgilar
        "gc_percent_all": 100 * (g + c) / total if total else 0.0,
        "gc_percent_acgt": 100 * (g + c) / acgt if acgt else 0.0,
    }


def main(paths):
    if not paths:
        sys.exit("Foydalanish: python3 gc_content.py <fasta> [<fasta> ...]")
    for p in paths:
        if not Path(p).is_file():
            print(f"Xato: fayl topilmadi: {p}", file=sys.stderr)
            continue
        for rec in SeqIO.parse(p, "fasta"):
            s = gc_stats(str(rec.seq))
            print(f"Fayl: {p}")
            print(f"  Yozuv: {rec.id}")
            print(f"  Uzunlik: {s['length']} bp "
                  f"(A={s['A']}, T={s['T']}, G={s['G']}, C={s['C']}, boshqa={s['other']})")
            print(f"  GC-tarkib (barcha belgilar): {s['gc_percent_all']:.2f}%")
            print(f"  GC-tarkib (faqat ACGT):      {s['gc_percent_acgt']:.2f}%")


if __name__ == "__main__":
    main(sys.argv[1:])
```

Ishga tushirish:

```bash
cd src/project1/script
python3 gc_content.py ../data/fasta1.fasta ../data/fasta2.fasta | tee ../results/3.3.txt
```

**Skriptning ishlash mantig'i (P2P uchun):**

1. `SeqIO.parse(..., "fasta")` fayldan har bir yozuvni o'qiydi (sarlavha va ketma-ketlikni ajratadi; qator bo'linishlarini o'zi birlashtiradi).
2. `upper()` — kichik harflarni bir xil ko'rinishga keltiradi.
3. `count()` bilan G, C (va A, T) sanaladi.
4. Ikki xil foiz chiqariladi: barcha belgilarga nisbatan va faqat ACGT'ga nisbatan — N belgilari bo'lsa farqi ko'rinadi.
5. Murakkablik: O(n) vaqt, n — ketma-ketlik uzunligi.

### 4.3. Etalon GC va P2P javobi (3.4-band)

**Savol:** «pUC19 vektorining etalon GC-tarkibi qancha? Olingan qiymatlar bilan solishtir. Fikringcha, farq bormi va nima uchun?»

**Javob asosi (raqamlarni o'zingnikiga almashtir):**

> **Etalon:** loyiha materiallarida pUC19 uchun ~51% keltirilgan; adabiyotda va ketma-ketlik bazalarida odatda 51–52% atrofida beriladi (aniq qiymat ketma-ketlik versiyasiga bog'liq). Mening natijalarim: `fasta1` — X.XX%, `fasta2` — Y.YY%.
>
> **Farq bormi?** *(skript natijasiga qarab javob ber)*
>
> Agar qiymatlar 1% dan kam farq qilsa: ular etalonga mos va ketma-ketliklar asosan bir xil vektor. Kichik farq bir necha nukleotid o'zgarishidan (variant, SNP) bo'lishi mumkin.
>
> Agar farq sezilarli bo'lsa, mumkin bo'lgan sabablar:
> 1. **Ketma-ketliklar uzunligi turlicha** — `fasta2` pUC19'ning o'zgartirilgan varianti (masalan, MCS'ga insert kiritilgan yoki kassetta qo'shilgan) bo'lishi mumkin; qo'shilgan DNK o'z GC-tarkibiga ega bo'lgani uchun umumiy foizni siljitadi.
> 2. **Boshqa tur/klon:** «pUC19» deb nomlangan yozuvlar ham har xil submissionlar bo'lib, ba'zilarida qo'shimcha elementlar yoki ketma-ketlik xatolari bo'ladi.
> 3. **Aniqlanmagan belgilar (N, IUPAC):** agar maxrajga qo'shilsa, foiz pasayadi. Shuning uchun ikkala variantni hisobladim.
> 4. **Hisoblash usuli:** barcha belgilar yoki faqat ACGT bo'yicha, registr (katta/kichik harf) bilan ishlash.
> 5. **Boshlanish nuqtasi va yo'nalish** (halqasimon molekulani qayerdan «kesib» yozish) GC foiziga ta'sir qilmaydi, chunki bu faqat ulushni hisoblaydi — buni ham eslatish yaxshi.

> **Eslatma:** `fasta1`/`fasta2` aniq nima ekanini (uzunligi, sarlavhasi) skript chiqishidagi `Yozuv` va `Uzunlik` qatorlaridan ko'r va NCBI sahifasidagi tavsifni o'qi — javobingdagi farq sababini shu ma'lumotga tayanib yoz. Men bu yozuvlar ichidagi ketma-ketliklarni ko'rmaganim uchun aniq sababni ayta olmayman.

---

## 5. Topshirish (Git)

```bash
cd <repo>
git status                       # faqat src/ ichidagi kerakli fayllar bo'lsin
git add src/
git commit -m "Add pUC19 analysis: FastQC, Trimmomatic, GC content"
git push origin develop
```

**Katta fayllar haqida:** `.fastq` fayllar yuzlab MB bo'lishi mumkin. README «topshiriqlarda ko'rsatilgan fayllardan boshqa fayl bo'lmasin» deydi, shuning uchun ma'lumotni repoga qo'shish-qo'shmaslikni P2P qoidalari/chat orqali aniqlashtir. Aniq bo'lmasa, `results/` dagi hisobotlar va kichik fayllarni commit qil.

---

## 6. Yakuniy tekshiruv ro'yxati

**Struktura**
- [ ] `src/project1/{data,results,script}` mavjud
- [ ] `develop` branch'ga push qilingan
- [ ] Ortiqcha fayllar yo'q

**1-topshiriq**
- [ ] `.fastq` `data/`da
- [ ] FastQC natijasi (`1.1.txt`) `results/`da
- [ ] Qizil/to'q sariq modullar belgilangan
- [ ] Adapterlar nomlari yozilgan

**2-topshiriq**
- [ ] Eng mos adapter fayli `data/`da
- [ ] Tozalangan fayl va FastQC natijalari `results/`da (2.2, 2.3, 2.5)
- [ ] Standart parametrlar bilan ikkinchi tozalash (2.7)
- [ ] Ikkala savolga (2.6, 2.8) javob tayyor

**3-topshiriq**
- [ ] `fasta1.fasta`, `fasta2.fasta` `data/`da
- [ ] `gc_content.py` `script/`da va ishlaydi
- [ ] Etalon GC bilan solishtirish va farq sababi tayyor (3.4)

**P2P**
- [ ] Har bir parametrni, har bir skript qatorini tushuntira olaman
