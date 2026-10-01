# pUC19 loyihasi uchun nazariy qo'llanma

Bu qo'llanma «pUC19 plazmidasini tahlil qilish» loyihasini bajarish uchun kerak bo'ladigan barcha tushuncha va ko'nikmalarni qadam-baqadam tushuntiradi. Texnik terminlar xalqaro (inglizcha) shaklda saqlangan.

## Mundarija

0. [Nimalarni o'rganishing kerak: xarita](#0-nimalarni-organishing-kerak-xarita)
1. [Biologiya asoslari: DNK va plazmida](#1-biologiya-asoslari-dnk-va-plazmida)
2. [Sekvenirlash va FASTQ format](#2-sekvenirlash-va-fastq-format)
3. [Phred Score (Q-score)](#3-phred-score-q-score)
4. [FastQC: hisobotni o'qish](#4-fastqc-hisobotni-oqish)
5. [Trimmomatic: tozalash](#5-trimmomatic-tozalash)
6. [GC-tarkib: nazariya, pseudocode, murakkablik, kod](#6-gc-tarkib)
7. [Amaliy ko'nikmalar: Linux, Git, Python, Java](#7-amaliy-konikmalar)
8. [Keyingi bosqichlar: alignment, SAM/BAM, VCF](#8-keyingi-bosqichlar)
9. [P2P'ga tayyorgarlik ro'yxati](#9-p2pga-tayyorgarlik)

---

## 0. Nimalarni o'rganishing kerak: xarita

| # | Mavzu | Nima uchun kerak | Chuqurlik |
|---|-------|------------------|-----------|
| 1 | DNK, nukleotid, plazmida, pUC19 tuzilishi | Loyiha kontekstini tushunish | O'rta |
| 2 | Sekvenirlash (Illumina), FASTQ, FASTA formatlari | Fayllarni o'qiy olish | Yuqori |
| 3 | Phred Score | Sifat nazoratining asosi | Yuqori |
| 4 | FastQC modullari (PASS/WARN/FAIL) | 1-topshiriq | Yuqori |
| 5 | Adapter, ILLUMINACLIP, SLIDINGWINDOW, TRAILING, MINLEN | 2-topshiriq | Yuqori |
| 6 | GC-tarkib va uni hisoblash | 3-topshiriq | Yuqori |
| 7 | Linux terminali (`cd`, `ls`, `head`, `wc`, `gunzip`) | Hamma narsa terminalda bajariladi | O'rta |
| 8 | Java (faqat ishga tushirish) | Trimmomatic — `.jar` fayl | Past |
| 9 | Python + Biopython (`SeqIO`) | 3-topshiriq skripti | O'rta |
| 10 | Git: `clone`, `checkout -b develop`, `add`, `commit`, `push` | Loyihani topshirish | O'rta |
| 11 | SRA ma'lumotlarini yuklash (SRA Toolkit yoki ENA) | 1-topshiriq, 2-band | O'rta |
| 12 | P2P'da o'z natijangni og'zaki tushuntirish | Himoya | Yuqori |

**Tavsiya etilgan tartib:** 1 → 2 → 3 → 4 → 5 → 6 → 7–11 (amaliyot bilan parallel) → 12.

---

## 1. Biologiya asoslari: DNK va plazmida

### 1.1. DNK va nukleotidlar

DNK — to'rt xil «harf»dan iborat uzun ketma-ketlik. Bu harflar nukleotidlar deyiladi:

| Harf | Nomi | Turi |
|------|------|------|
| **A** | Adenin | Purin |
| **T** | Timin | Pirimidin |
| **G** | Guanin | Purin |
| **C** | Sitozin | Pirimidin |

DNK ikki zanjirdan iborat (double helix) va ular **komplementarlik** qoidasi bilan bog'langan: **A–T** (2 ta vodorod bog'i) va **G–C** (3 ta vodorod bog'i). G–C juftligi 3 ta bog'ga ega bo'lgani uchun A–T'dan mustahkamroq — bu GC-tarkib ahamiyatining sabablaridan biri (6-bo'limga qarang).

**bp (base pair)** — nukleotid juftligi, DNK uzunligining o'lchov birligi. Masalan, «2686 bp» — 2686 ta juft.

### 1.2. Plazmida

**Plazmida** — bakteriya ichidagi xromosomadan alohida joylashgan, o'z-o'zidan ko'payadigan kichik halqasimon DNK molekulasi. Biologlar plazmidalarni «tashuvchi» (vektor) sifatida ishlatadi: kerakli DNK fragmentini plazmidaga joylashtirib, bakteriyaga kiritadi, bakteriya esa uni ko'paytiradi.

### 1.3. pUC19 tuzilishi

pUC19 — eng mashhur klonlash vektorlaridan biri. Uzunligi taxminan **2686 bp**. Asosiy elementlari:

| Element | Vazifasi |
|---------|----------|
| **ori** (origin of replication) | Plazmida ko'payishi boshlanadigan joy. Mutatsiya bo'lsa, kopiyalar soni (copy number) pasayishi mumkin |
| **AmpR** (β-laktamaza geni) | Ampitsillin antibiotigiga chidamlilik. Faqat plazmida kirgan bakteriyalar antibiotikli muhitda yashaydi |
| **lacZα** | β-galaktozidaza fermentining bir qismini kodlaydi (ko'k-oq selektsiya uchun) |
| **MCS** (multiple cloning site) | Ko'plab restriksiya fermentlari kesadigan qisqa hudud — kerakli DNK (insert) shu yerga joylashtiriladi |

### 1.4. Ko'k-oq selektsiya (blue-white screening)

1. MCS **lacZα** genining ichida joylashgan.
2. Insert **yo'q** bo'lsa → lacZα ishlaydi → X-gal substrati ko'k rangga aylanadi → **ko'k koloniya**.
3. Insert **bor** bo'lsa → lacZα buziladi → rang hosil bo'lmaydi → **oq koloniya**.

Demak, oq koloniya = «ehtimol insert bor». Lekin bu faqat dastlabki saralash: insert to'g'rimi, to'liqmi, mutatsiya yo'qmi — buni faqat **sekvenirlash** aniq ko'rsatadi. README'dagi butun «bir kun» shu savolga javob izlashdir.

---

## 2. Sekvenirlash va FASTQ format

### 2.1. Sekvenirlash nima?

Sekvenator (loyihada — Illumina turidagi) DNK'ni o'qiydi, lekin butun molekulani bir martada emas. Jarayon taxminan quyidagicha:

1. DNK minglab nusxada olinadi.
2. Mayda parchalarga bo'linadi (fragmentation).
3. Parchalarning uchlariga sun'iy **adapter** ketma-ketliklari «tikiladi» (ular sekvenator uchun «tutqich» vazifasini bajaradi).
4. Har bir parcha o'qiladi → natija **read** (o'qilma) deyiladi. Read — odatda 50–300 nukleotidlik qisqa ketma-ketlik.

**Single-end (SE)** — parcha faqat bir uchidan o'qiladi (1 ta fayl). **Paired-end (PE)** — ikkala uchidan o'qiladi (2 ta fayl: `_1.fastq` va `_2.fastq`). Trimmomatic'da `SE` yoki `PE` rejimini ma'lumotingga qarab tanlaysan — bu juda muhim.

**Coverage (qoplama chuqurligi)** — genomning har bir pozitsiyasi o'rtacha nechta read bilan «yopilgani».

### 2.2. FASTQ format

FASTQ faylida har bir read **aniq 4 qatorda** yoziladi:

```
@SRR1553610.1 HWI-ST:1:1101:1234:2000 length=36    ← 1) Identifikator (@ bilan boshlanadi)
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAA              ← 2) Nukleotid ketma-ketlik
+                                                  ← 3) Ajratuvchi (+)
IIIIIIIIIHIIIIIGIIIIFDDDD@@@@?5555;;              ← 4) Har bir nukleotid uchun sifat (ASCII)
```

**Qoida:** 2-qator va 4-qator uzunligi doim teng. 4-qatordagi har bir belgi 2-qatordagi mos nukleotidning sifatini kodlaydi.

Faylda read'lar soni = (qatorlar soni) / 4. Buni terminalda tekshirish:

```bash
echo $(( $(wc -l < data/reads.fastq) / 4 ))
```

### 2.3. FASTA format

FASTA — bitta yoki bir nechta ketma-ketlikni saqlaydigan oddiy format (sifat ma'lumoti **yo'q**). 3-topshiriqda aynan shu formatda ishlaysan:

```
>L09137.2 Cloning vector pUC19
TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCA
CAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTG
...
```

`>` bilan boshlanuvchi qator — sarlavha (header), keyingi qatorlar — ketma-ketlik (odatda 60–80 belgidan qatorlarga bo'lingan).

| | FASTQ | FASTA |
|---|---|---|
| Nima saqlaydi | Read + sifat | Ketma-ketlik |
| Qatorlar/yozuv | 4 | 1 sarlavha + N qator |
| Odatda | Xom sekvenirlash ma'lumoti | Referens, tayyor ketma-ketlik |

---

## 3. Phred Score (Q-score)

### 3.1. Formula

Sekvenator har bir nukleotidni o'qiganda xato qilishi mumkin. Phred sifat balli xato ehtimolini logarifmik shkalada ifodalaydi:

```
Q = -10 · log10(P)        ⇔        P = 10^(-Q/10)
```

bu yerda **P** — shu nukleotid noto'g'ri aniqlangan bo'lish ehtimoli.

| Q-score | Xato ehtimoli P | Aniqlik | Ma'nosi |
|---------|-----------------|---------|---------|
| 10 | 1/10 = 10% | 90% | Yomon |
| 20 | 1/100 = 1% | 99% | Qoniqarli |
| 30 | 1/1000 = 0.1% | 99.9% | Yaxshi |
| 40 | 1/10 000 = 0.01% | 99.99% | A'lo |

Qoida: Q har **10** ga oshganda xato **10 baravar** kamayadi.

### 3.2. ASCII kodlash (Phred+33)

4-qatordagi belgi sifatni ASCII kodi orqali bildiradi:

```
belgi → ASCII kodi → Q = kod − 33
```

Masalan: `I` → ASCII 73 → Q = 73 − 33 = **40**. `5` → ASCII 53 → Q = **20**. `#` → ASCII 35 → Q = **2** (juda yomon).

Bu «Phred+33» (Sanger / zamonaviy Illumina 1.8+) kodlashi. Eski Illumina 1.3–1.5 «Phred+64» ishlatgan. Trimmomatic'da `-phred33` yoki `-phred64` parametri shu uchun kerak. Zamonaviy ma'lumotlarda deyarli har doim `-phred33`.

### 3.3. Python'da hisoblash

```python
def phred_from_char(ch: str, offset: int = 33) -> int:
    return ord(ch) - offset

print(phred_from_char("I"))   # 40
print(phred_from_char("5"))   # 20
```

---

## 4. FastQC: hisobotni o'qish

FastQC FASTQ faylni tahlil qilib, HTML hisobot va `.zip` arxiv (ichida `summary.txt`, `fastqc_data.txt`) yaratadi. Har bir modul **PASS / WARN / FAIL** status oladi.

> **Muhim:** FAIL = «albatta muammo» degani emas. FastQC umumiy (genom) sekvenirlash uchun mo'ljallangan, ba'zi biologik holatlarda (masalan, plazmida) ba'zi FAIL'lar kutilgan bo'lishi mumkin. Sen buni **tushuntira olishing** kerak.

### 4.1. Modullar jadvali

| Modul | Nimani ko'rsatadi | Yaxshi holat | Muammo belgisi va sabablari |
|-------|-------------------|--------------|------------------------------|
| **Basic Statistics** | Read'lar soni, uzunligi, GC%, kodlash | Odatda PASS | Kamdan-kam FAIL beradi |
| **Per base sequence quality** | Har bir pozitsiyadagi sifat (box-plot) | Hammasi yashil zonada (Q>28) | Oxirga qarab sifat tushishi — sekvenator «charchashi» |
| **Per tile sequence quality** | Flow cell'ning har bir qismi (tile) sifati | Bir tekis ko'k | Qizil dog'lar — texnik muammo (pufak, dog') |
| **Per sequence quality scores** | Read'larning o'rtacha sifat taqsimoti | Yuqori Q'da bitta cho'qqi | Past Q'da alohida cho'qqi — read'larning bir qismi yomon |
| **Per base sequence content** | Har pozitsiyada A/T/G/C ulushi | 4 ta chiziq deyarli parallel | Boshida notekislik — kutubxona tayyorlash artefakti yoki adapterlar |
| **Per sequence GC content** | Read'lar GC% taqsimoti vs nazariy normal taqsimot | Silliq qo'ng'iroq shakli | Keng/ikki cho'qqili — ifloslanish. Plazmidada kichik genom tufayli «tishli» bo'lishi mumkin |
| **Per base N content** | Aniqlanmagan (N) nukleotidlar ulushi | ~0% | N ko'p — sekvenator o'qiy olmagan |
| **Sequence length distribution** | Read uzunliklari taqsimoti | Bitta cho'qqi (barcha read bir xil uzunlikda) | Trimmomatic'dan keyin turli uzunliklar — **kutilgan holat** |
| **Sequence duplication levels** | Takrorlangan read'lar ulushi | Past | Yuqori — PCR duplikatlari yoki kichik genomning juda yuqori coverage'i |
| **Overrepresented sequences** | Odatdagidan ko'p uchraydigan ketma-ketliklar (>0.1%) | Ro'yxat bo'sh | Adapterlar, ifloslanish. FastQC ba'zan «Possible Source» ustunida adapter nomini ko'rsatadi |
| **Adapter content** | Adapter ketma-ketligi bor read'lar ulushi (pozitsiya bo'yicha) | Chiziq 0% da | Chiziq 5–10% dan oshsa — WARN/FAIL; adapter nomi grafikda legendada ko'rinadi |

### 4.2. 1-topshiriq uchun amaliy maslahat

- **Adapter nomini** topish: *Adapter Content* grafigi legendasi (masalan, «Illumina Universal Adapter», «Nextera Transposase Sequence») va *Overrepresented sequences* jadvalining «Possible Source» ustuni.
- **Qizil va to'q sariq** statuslarni aniq `summary.txt` faylidan olish mumkin — har qator `PASS/WARN/FAIL  Modul nomi  fayl nomi` ko'rinishida:

```bash
grep -E "^(FAIL|WARN)" results/*_fastqc/summary.txt
```

### 4.3. Ishga tushirish

```bash
fastqc data/reads.fastq -o results/
```

`-o` — natija papkasi (u oldindan mavjud bo'lishi kerak).

---

## 5. Trimmomatic: tozalash

### 5.1. Nima uchun tozalash kerak?

- **Adapterlar** — DNK'ga emas, sekvenirlash texnologiyasiga tegishli sun'iy ketma-ketliklar. Read oxirida paydo bo'lsa, keyingi alignment/assembly'ni buzadi.
- **Past sifatli uchlar** — read oxirida Q pasayadi; bu yerdagi «mutatsiya»lar aslida sekvenator xatosi bo'lishi mumkin va soxta variantlar (false positive SNP) keltirib chiqaradi.
- **Juda qisqa read'lar** — noaniq joylashadi (ko'p joyga mos tushadi), ma'lumot sifatini pasaytiradi.

### 5.2. Ishga tushirish sintaksisi

**Single-end:**

```bash
java -jar trimmomatic-0.39.jar SE -phred33 \
    input.fastq output.fastq \
    ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
```

**Paired-end:**

```bash
java -jar trimmomatic-0.39.jar PE -phred33 \
    in_1.fastq in_2.fastq \
    out_1P.fastq out_1U.fastq out_2P.fastq out_2U.fastq \
    ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:2:True \
    SLIDINGWINDOW:4:15 MINLEN:36
```

PE rejimida 4 ta chiqish fayli: `P` (paired — juft saqlanganlar) va `U` (unpaired — juftidan ayrilganlar).

> **Muhim:** qadamlar **yozilgan tartibda** bajariladi. ILLUMINACLIP odatda **birinchi** turishi kerak (adapterni kesmay turib sifat bo'yicha kesish adapterni topishni qiyinlashtiradi).

### 5.3. Funksiyalar

| Funksiya | Parametrlar | Nima qiladi |
|----------|-------------|-------------|
| **ILLUMINACLIP** | `fasta:seedMismatches:palindromeThreshold:simpleClipThreshold[:minAdapterLength:keepBothReads]` | Adapterlarni topib kesadi |
| **LEADING** | `quality` | Read **boshidan** sifati ostonadan past nukleotidlarni kesadi |
| **TRAILING** | `quality` | Read **oxiridan** sifati ostonadan past nukleotidlarni kesadi |
| **SLIDINGWINDOW** | `windowSize:requiredQuality` | Oyna bo'ylab yuradi; oynadagi o'rtacha sifat ostonadan tushsa, shu joydan keyingisini kesadi |
| **MINLEN** | `length` | Kesishdan keyin shu uzunlikdan qisqa read'ni butunlay tashlaydi |
| **CROP** | `length` | Read'ni berilgan uzunlikka (oxiridan) qisqartiradi |
| **HEADCROP** | `length` | Boshidan N ta nukleotidni kesadi |
| **AVGQUAL** | `quality` | O'rtacha sifati past read'ni tashlaydi |
| **TOPHRED33 / TOPHRED64** | — | Sifat kodlashini o'zgartiradi |

### 5.4. ILLUMINACLIP parametrlarini tushunish

README'da «2 mismatches, 30, 10» deb yozilgan. Aniqroq aytganda, Trimmomatic manualiga ko'ra:

- **seedMismatches (2)** — adapterning qisqa «urug'i» (seed, ~16 nt) read bilan nechta nomuvofiqlik bilan mos kelishi mumkin. Katta bo'lsa — ko'proq adapter topiladi, lekin soxta topilmalar ham ko'payadi.
- **palindromeClipThreshold (30)** — faqat **PE** rejimida: ikki read adapter orqali «palindrom» bo'lib mos kelganda qo'llaniladi. Yuqori qiymat — qat'iyroq.
- **simpleClipThreshold (10)** — adapter read bilan «oddiy» (bitta read bo'yicha) mos kelishining ostonasi. **Faqat shu SE rejimida ham ishlaydi.**
- **minAdapterLength** va **keepBothReads** — faqat PE palindrom rejimi uchun.

**Eslatma:** parametrlar tartibi `seedMismatches : palindromeThreshold : simpleClipThreshold`. Agar sen README'dagi (2, 30, 10) ni kiritsang, bu `2:30:10` bo'ladi. SE rejimida palindrom ostonasi ishlatilmaydi, lekin sintaksis uchun baribir yoziladi. O'z javobingda shuni to'g'ri tushuntirishing kerak — bu P2P'da so'raladigan nozik joy.

### 5.5. Adapter fayllarini tanlash

Trimmomatic bilan keladigan `adapters/` papkasida standart fayllar bor:

| Fayl | Qachon |
|------|--------|
| `TruSeq3-SE.fa` | Illumina HiSeq/MiSeq, single-end |
| `TruSeq3-PE.fa` / `TruSeq3-PE-2.fa` | HiSeq/MiSeq, paired-end |
| `TruSeq2-SE.fa` / `TruSeq2-PE.fa` | Eski Illumina GAII |
| `NexteraPE-PE.fa` | Nextera kutubxonalari |

Qaysi birini tanlashni **FastQC'dagi adapter nomi** va **sekvenirlash tavsifi** (NCBI SRA sahifasidagi Platform, Library Layout, Instrument) hal qiladi. Adapter yo'qolmasa — boshqa faylni sinab ko'r (2-topshiriq, 4-band aynan shuni so'raydi).

### 5.6. Nima uchun read uzunligi taqsimoti o'zgaradi?

Xom ma'lumotda barcha read'lar odatda **bir xil uzunlikda** (masalan, 100 nt) — shuning uchun *Sequence length distribution* — bitta cho'qqi. Trimmomatic har bir read'ni **alohida** kesadi: biriga adapter bor (qisqaradi), birining oxiri sifatsiz (qisqaradi), boshqasi umuman o'zgarmaydi. Natijada uzunliklar turlicha bo'lib qoladi → grafik «keng» bo'ladi va FastQC WARN berishi mumkin. Bu **normal va kutilgan** natija, xato emas. (Bu 2-topshiriq 6-bandi javobining asosi.)

---

## 6. GC-tarkib

### 6.1. Nazariya

**GC-tarkib** (GC content) — ketma-ketlikdagi G va C nukleotidlar ulushi:

```
GC% = (G + C) / (A + T + G + C) × 100
```

**Nima uchun muhim?**

1. **Termostabillik:** G–C juftligida 3 ta, A–T da 2 ta vodorod bog'i bor → GC ko'p DNK yuqori haroratda erishga (melting) chidamli.
2. **Tur belgisi:** har bir organizmning genomi xarakterli GC-tarkibga ega (odam ~41%, *E. coli* ~51%). Ifloslanishni aniqlashda foydali.
3. **Sekvenirlash sifati:** juda yuqori yoki past GC hududlar sekvenatorda notekis qoplanadi (GC bias).
4. **Tezkor tekshiruv:** README'dagi kabi — kutilgan GC bilan solishtirib, namuna to'g'riligini taxminan baholash.

**pUC19 uchun** README etalon sifatida ~51% ni beradi. (Aniq qiymat qaysi ketma-ketlik versiyasi va qanday hisoblanishiga qarab 51–52% atrofida bo'lishi mumkin; buni o'zing hisoblab ko'r.)

### 6.2. Hisoblashdagi nozik joylar

- **Registr:** `g` va `G` bir xil hisoblanishi kerak → `upper()`.
- **N va boshqa IUPAC belgilar:** `N` — «istalgan nukleotid». Uni maxrajga qo'shish yoki qo'shmaslik natijani o'zgartiradi. Ikkala variantni hisoblab, farqni ko'rsatish yaxshi.
- **Ko'p yozuvli FASTA:** fayl ichida bir nechta ketma-ketlik bo'lishi mumkin — har birini alohida hisoblash kerak.
- **Sarlavha va yangi qator:** `>` qatorlar va `\n` hisobga kirmasligi kerak.

### 6.3. Pseudocode

```
FUNCTION GCContent(sequence):
    seq ← UPPERCASE(sequence)
    gc ← 0
    total ← 0
    FOR EACH nucleotide IN seq:
        total ← total + 1
        IF nucleotide = 'G' OR nucleotide = 'C':
            gc ← gc + 1
    IF total = 0:
        RETURN 0
    RETURN (gc / total) × 100
```

Fayl bilan ishlovchi to'liq variant:

```
FUNCTION GCContentFromFasta(path):
    OPEN file path
    FOR EACH record IN file:                 // '>' sarlavha — yangi yozuv
        sequence ← JOIN(all non-header lines of record)
        PRINT record.id, GCContent(sequence)
```

### 6.4. Murakkablik

- **Vaqt:** O(n), n — ketma-ketlik uzunligi. Har bir nukleotid bir marta ko'riladi (pUC19 uchun n ≈ 2686, deyarli bir lahza).
- **Xotira:** satrma-satr o'qilsa O(1) qo'shimcha (faqat sanagichlar). Butun ketma-ketlik xotirada saqlansa — O(n). Katta genomlar (milliardlab nt) uchun satrma-satr o'qish zarur; pUC19 uchun farqi yo'q.
- `str.count("G") + str.count("C")` ham O(n) (ikki o'tish), amalda C darajasida bajarilgani uchun Python'da tez.

### 6.5. Implementatsiyalar

#### Python (Biopython)

```python
from Bio import SeqIO

def gc_percent(sequence: str) -> float:
    seq = sequence.upper()
    total = len(seq)
    if total == 0:
        return 0.0
    return 100 * (seq.count("G") + seq.count("C")) / total

for record in SeqIO.parse("data/fasta1.fasta", "fasta"):
    print(record.id, len(record.seq), f"{gc_percent(str(record.seq)):.2f}%")
```

Biopython'da tayyor funksiya ham bor (`Bio.SeqUtils.gc_fraction`), lekin avval qo'lda yozib tushunib ol — P2P'da «qanday ishlaydi?» deb so'rashadi.

#### C++

```cpp
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 2) { std::cerr << "Foydalanish: " << argv[0] << " <fasta>\n"; return 1; }
    std::ifstream in(argv[1]);
    if (!in) { std::cerr << "Fayl ochilmadi\n"; return 1; }
    std::string line;
    long long gc = 0, acgt = 0, total = 0;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') continue;      // header
        for (unsigned char ch : line) {
            if (std::isspace(ch)) continue;
            char u = std::toupper(ch);
            ++total;
            if (u == 'G' || u == 'C') { ++gc; ++acgt; }
            else if (u == 'A' || u == 'T') ++acgt;
        }
    }
    std::cout << std::fixed << std::setprecision(2)
              << "Uzunlik: " << total << " bp\n"
              << "GC (barcha): " << (total ? 100.0 * gc / total : 0.0) << "%\n"
              << "GC (ACGT):   " << (acgt ? 100.0 * gc / acgt : 0.0) << "%\n";
}
```

Kompilyatsiya: `g++ -O2 -o gc gc.cpp && ./gc data/fasta1.fasta`

#### C

```c
#include <ctype.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Foydalanish: %s <fasta>\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("fopen"); return 1; }
    char buf[65536];
    long long gc = 0, acgt = 0, total = 0;
    while (fgets(buf, sizeof buf, f)) {
        if (buf[0] == '>') continue;                 /* header */
        for (char *p = buf; *p; ++p) {
            if (isspace((unsigned char)*p)) continue;
            int u = toupper((unsigned char)*p);
            ++total;
            if (u == 'G' || u == 'C') { ++gc; ++acgt; }
            else if (u == 'A' || u == 'T') ++acgt;
        }
    }
    fclose(f);
    printf("Uzunlik: %lld bp\n", total);
    printf("GC (barcha): %.2f%%\n", total ? 100.0 * gc / total : 0.0);
    printf("GC (ACGT):   %.2f%%\n", acgt ? 100.0 * gc / acgt : 0.0);
    return 0;
}
```

Kompilyatsiya: `gcc -O2 -o gc gc.c && ./gc data/fasta1.fasta`

> C va C++ variantlari ixtiyoriy (loyiha Python'ni talab qiladi), lekin C/C++ o'rganayotgan bo'lsang, mashq sifatida foydali. Uchala kod ham sinab ko'rilgan.

#### Misol bo'yicha qo'lda tekshiruv

`GGCCAATT` → G=2, C=2, A=2, T=2 → GC = 4/8 = **50%**. Skriptingni avval shunday kichik misolda sinab ko'r.

---

## 7. Amaliy ko'nikmalar

### 7.1. Linux terminali (minimum to'plam)

```bash
mkdir -p project1/{results,data,script}   # papkalar strukturasi (bir buyruqda)
cd project1                                # papkaga kirish
ls -lh data/                               # fayllar va o'lchamlari
head -n 8 data/reads.fastq                 # birinchi 2 ta read (8 qator)
wc -l data/reads.fastq                     # qatorlar soni
gunzip data/reads.fastq.gz                 # .gz arxivni ochish
grep -c "^>" data/fasta1.fasta             # FASTA'dagi yozuvlar soni
cp /yo'l/TruSeq3-SE.fa data/               # nusxalash
```

### 7.2. Java

Trimmomatic — `.jar` fayl, Java kerak:

```bash
java -version          # o'rnatilganini tekshirish
# Ubuntu/Debian: sudo apt install default-jre
```

Yoki conda orqali hammasini bir yo'la o'rnatish:

```bash
conda install -c bioconda fastqc trimmomatic sra-tools
```

(conda'da `trimmomatic` buyrug'i to'g'ridan-to'g'ri ishlaydi, `java -jar` shart emas.)

### 7.3. SRA'dan ma'lumot yuklash

README sahifadan qo'lda yuklashni aytadi (NCBI → Run → FASTQ). Buyruq satri varianti:

```bash
prefetch SRRxxxxxxx
fasterq-dump SRRxxxxxxx -O data/
```

`SRRxxxxxxx` o'rniga NCBI sahifasidagi **Run** identifikatorini qo'y. Ma'lumot paired-end bo'lsa, ikkita fayl (`_1`, `_2`) hosil bo'ladi — bu SE/PE tanlovingga ta'sir qiladi, shuning uchun NCBI sahifasidagi **Layout** maydonini albatta o'qi.

### 7.4. Git (loyiha talabi)

```bash
git clone <GitLab-repo-URL>
cd <repo>
git checkout -b develop            # develop branch yaratish
# fayllarni src/ ichida yarat
git add src/
git commit -m "Add pUC19 analysis"
git push origin develop            # aynan develop'ni push qil
```

**Eslatma:** README talabi — `src/` papkada topshiriqlarda ko'rsatilgandan **boshqa fayl bo'lmasligi** kerak. Katta `.fastq` fayllarni commit qilma (o'lchami katta, `.gitignore`ga qo'sh) — avval P2P qoidalarida shu haqda ko'rsatma bor-yo'qligini tekshir.

### 7.5. Python muhiti

```bash
python3 --version             # 3.8+
pip install biopython
python3 -c "import Bio; print(Bio.__version__)"
```

---

## 8. Keyingi bosqichlar

README hikoyasida bu loyihada **bajarilmaydigan**, lekin keyingi loyihalarda uchraydigan bosqichlar ham tilga olingan. Ular nima ekanini bilib qo'yish P2P'da foydali:

| Bosqich | Vosita | Nima qiladi |
|---------|--------|-------------|
| **Alignment** (hizalash) | BWA, Bowtie2 | Har bir read'ni referens ketma-ketlik bo'yicha joylashtiradi |
| **SAM/BAM** | samtools | Alignment natijalari (BAM — SAM'ning siqilgan binar shakli) |
| **Variant calling** | GATK | Referensdan farqlarni (SNP, Indel) topadi |
| **VCF** | — | Topilgan variantlar hisobot fayli |

- **SNP** — bitta nukleotid almashinuvi (A→T).
- **Indel** — kichik insertsiya yoki deletsiya.
- Haqiqiy mutatsiyani sekvenator xatosidan ajratish uchun **Q-score** va **coverage depth** qaraladi: bir nechta read bir xil farqni ko'rsatsa — haqiqiy bo'lish ehtimoli yuqori.

Mantiqiy zanjir: **FASTQ → FastQC → Trimmomatic → FastQC → alignment → variant calling → xulosa**. Sen hozir ushbu zanjirning birinchi uchdan bir qismini bajarasan.

---

## 9. P2P'ga tayyorgarlik

O'zingdan so'rab ko'r — har biriga 1–2 jumlada javob bera olasanmi?

- [ ] FASTQ'ning 4 qatori nimani anglatadi?
- [ ] Q=20 va Q=30 farqi nima? Q=30 uchun xato ehtimoli qancha?
- [ ] Adapter nima va nima uchun read'larda paydo bo'ladi?
- [ ] FastQC'da FAIL har doim yomonmi? Misol keltir.
- [ ] Qaysi adapterlar topildi va ularni qayerdan bilding (grafik / Overrepresented jadvali)?
- [ ] Nima uchun aynan shu adapter faylini tanlading?
- [ ] ILLUMINACLIP'dagi har bir raqam nimani bildiradi?
- [ ] SLIDINGWINDOW, TRAILING, MINLEN nima qiladi va nima uchun ularni shunday sozlading?
- [ ] Trimmomatic'dan keyin read uzunligi taqsimoti nima uchun o'zgardi?
- [ ] GC% formulasi, nima uchun muhim, skriptdagi har bir qator nima qiladi?
- [ ] Ikki fasta o'rtasida GC farqi bo'lsa — sabablari qanday bo'lishi mumkin?
- [ ] SE va PE farqi; sening ma'lumoting qaysi turda edi?

**Eng muhim qoida (README'dan):** boshqalarning yechimini ko'chirma. Har bir parametr va har bir qatorni tushuntira olishing kerak.
