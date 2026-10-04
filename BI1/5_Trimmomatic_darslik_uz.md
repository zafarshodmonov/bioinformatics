# Trimmomatic darsligi

> Manba: *Trimmomatic Manual V0.32* (Usadel Lab). Darslik o'zbek tilida, texnik terminlar inglizcha saqlangan. Oxirida sening **pUC19 loyihang** (`SRR3341836_1.fastq`) uchun amaliy qism bor.

## Mundarija

1. [Trimmomatic nima va nima uchun kerak](#1-trimmomatic-nima-va-nima-uchun-kerak)
2. [Asosiy tushunchalar](#2-asosiy-tushunchalar)
3. [Ishga tushirish: SE va PE rejimlari](#3-ishga-tushirish-se-va-pe-rejimlari)
4. [Kirish/chiqish fayllari](#4-kirishchiqish-fayllari)
5. [Bosqichlar (steps) batafsil](#5-bosqichlar-steps-batafsil)
6. [Adapter fayllari](#6-adapter-fayllari)
7. [Qo'llanmadagi misollar](#7-qollanmadagi-misollar)
8. [Algoritmlar: pseudocode, murakkablik, kod](#8-algoritmlar-pseudocode-murakkablik-kod)
9. [Sening loyihang uchun amaliy qism](#9-sening-loyihang-uchun-amaliy-qism)
10. [Mashqlar va javoblar](#10-mashqlar-va-javoblar)
11. [Qisqa ma'lumotnoma (cheat sheet)](#11-qisqa-malumotnoma)
12. [Versiyalar tarixi](#12-versiyalar-tarixi)

---

## 1. Trimmomatic nima va nima uchun kerak

**Trimmomatic** — Illumina (FASTQ) ma'lumotlarini kesish (trim) va qirqish (crop), shuningdek adapterlarni olib tashlash uchun tez, ko'p oqimli (multithreaded) buyruq satri dasturi.

Nima uchun kerak? Xom ma'lumotda ikki turdagi «axlat» bo'ladi:

1. **Adapterlar** — sekvenirlash uchun DNK uchlariga ulangan sun'iy ketma-ketliklar. Ular keyingi tahlilda (alignment, assembly) xalaqit beradi.
2. **Past sifatli asoslar** — read oxirida sekvenator xatolari ko'payadi, bular soxta mutatsiyalarga olib keladi.

Dastur ikki rejimda ishlaydi:

| Rejim | Kirish | Xususiyat |
|---|---|---|
| **SE** (single-end) | 1 ta fayl | Oddiy |
| **PE** (paired-end) | 2 ta fayl (forward, reverse) | Juft read'lar mosligini saqlaydi va adapter «read-through»ni topishda qo'shimcha ma'lumotdan foydalanadi |

Trimmomatic **phred+33** yoki **phred+64** kodlashli FASTQ bilan ishlaydi. `.gz` va `.bz2` bilan siqilgan fayllarni to'g'ridan-to'g'ri o'qiydi.

---

## 2. Asosiy tushunchalar

### 2.1. Adapter «read-through» (nima uchun adapter read oxirida paydo bo'ladi)

Qo'llanmaga ko'ra, adapter ifloslanishining **eng keng tarqalgan sababi** — read uzunligidan **qisqa** DNK parchasi sekvenirlanishi. Read boshida haqiqiy ma'lumot bo'ladi, parcha tugagach sekvenator adapter ichiga «o'qishda davom etadi» (read-through):

```
Parcha (DNK):   |=================|
Read (150 bp):  [===== DNK =====][--- adapter ---]
                 0            60 61              150
                                ↑ kesish nuqtasi
```

To'liq adapterni topish oson, qisman (qisqa) adapterni ishonchli topish esa qiyin. Bu 5.1-bo'limdagi sezgirlik/o'ziga xoslik kelishuvining (trade-off) sababi.

PE ma'lumotda read-through forward va reverse read'larda **bir xil pozitsiyada** sodir bo'ladi, va ularning adapterdan oldingi qismlari bir-birining **reverse-complement**i bo'ladi. Trimmomatic'ning «palindrome» rejimi aynan shuni ishlatadi.

### 2.2. Phred Score (qisqacha)

`Q = -10 · log10(P)`. FASTQ'da sifat ASCII belgi bilan yoziladi: `Q = ASCII − 33` (phred+33). Masalan: `I` → 40, `5` → 20, `+` → 10, `#` → 2.

Illumina'da **Q=2** maxsus ma'noga ega: u «low quality segment» (read oxiridagi ishonchsiz bo'lak) belgisi. Qo'llanma shu bo'laklarni TRAILING bilan olib tashlash mumkinligini, lekin buning o'rniga SLIDINGWINDOW yoki MAXINFO tavsiya etilishini aytadi.

### 2.3. Qayta ishlash tartibi (Processing Order)

Bosqichlar buyruq satrida **yozilgan tartibda** bajariladi. Qo'llanma tavsiyasi: adapter kesish (agar kerak bo'lsa) **imkon qadar erta** bajarilsin, chunki qisman moslik orqali adapterni to'g'ri aniqlash qiyinroq (sifat bo'yicha kesish adapterning bir qismini olib tashlashi mumkin).

MINLEN esa odatda **eng oxirida** turishi kerak, chunki u boshqa bosqichlar read'ni qisqartirgandan **keyin** uzunlikni tekshiradi.

---

## 3. Ishga tushirish: SE va PE rejimlari

### 3.1. Single-end (SE)

```bash
java -jar <trimmomatic.jar> SE [-threads <n>] [-phred33 | -phred64] [-trimlog <logFayl>] \
    <input> <output> <step 1> <step 2> ...
```

Bitta kirish, bitta chiqish fayli.

### 3.2. Paired-end (PE)

```bash
java -jar <trimmomatic.jar> PE [-threads <n>] [-phred33 | -phred64] [-trimlog <logFayl>] \
    [-basein <inputBase> | <input1> <input2>] \
    [-baseout <outputBase> | <paired1> <unpaired1> <paired2> <unpaired2>] \
    <step 1> ...
```

Ikki kirish va **to'rt** chiqish fayli:

| Chiqish | Ma'nosi |
|---|---|
| forward paired | ikkala read ham omon qolgan, forward |
| forward unpaired | forward omon qoldi, juftligi tashlandi |
| reverse paired | ikkala read ham omon qolgan, reverse |
| reverse unpaired | reverse omon qoldi, juftligi tashlandi |

### 3.3. Umumiy parametrlar

| Parametr | Ma'nosi |
|---|---|
| `-phred33` / `-phred64` | Sifat kodlash. **0.32 versiyadan** ko'rsatilmasa, avtomatik aniqlanadi (avval standart `-phred64` edi) |
| `-threads <n>` | Oqimlar soni. Ko'rsatilmasa, avtomatik tanlanadi |
| `-trimlog <fayl>` | Har bir read uchun kesish jurnali |

**Trimlog** quyidagilarni yozadi: read nomi, omon qolgan ketma-ketlik uzunligi, birinchi omon qolgan asos o'rni (boshidan nechta kesilgani), oxirgi omon qolgan asosning asl read'dagi o'rni, oxiridan nechta kesilgani.

> Trimlog fayli katta bo'ladi (har bir read uchun qator). Diskda joying kam bo'lsa, uni yoqma.

---

## 4. Kirish/chiqish fayllari

### 4.1. Siqish (gzip/bzip2)

Fayl nomiga `.gz` yoki `.bz2` qo'shsang, Trimmomatic kirishni siqilgan deb o'qiydi, chiqishni esa o'zi siqib yozadi:

```bash
java -jar trimmomatic.jar SE input.fastq.gz output.fastq.gz ...
```

> **Amaliy foyda:** diskingda joy kam bo'lsa (VirtualBox'dagi 25 GB), `.gz` chiqish fayllari bir necha marta kichik bo'ladi. FastQC ham `.gz` o'qiy oladi.

### 4.2. PE uchun nomlash

**Kirish:** ikkita fayl nomini aniq yoz yoki `-basein` bilan faqat forward faylni ko'rsat; reverse fayl nomlash andozalaridan topiladi:

- `Sample_Name_R1_001.fq.gz` → `Sample_Name_R2_001.fq.gz`
- `Sample_Name.f.fastq` → `Sample_Name.r.fastq`
- `Sample_Name.1.sequence.txt` → `Sample_Name.2.sequence.txt`

**Chiqish:** to'rt nomni aniq yoz yoki `-baseout` bilan bitta nom ber. `-baseout mySampleFiltered.fq.gz` uchun:

| Fayl | Mazmuni |
|---|---|
| `mySampleFiltered_1P.fq.gz` | paired forward |
| `mySampleFiltered_1U.fq.gz` | unpaired forward |
| `mySampleFiltered_2P.fq.gz` | paired reverse |
| `mySampleFiltered_2U.fq.gz` | unpaired reverse |

---

## 5. Bosqichlar (steps) batafsil

Ko'p bosqichlar sozlamalarni **ikki nuqta (`:`)** bilan ajratib oladi.

### 5.1. ILLUMINACLIP

Illumina adapterlari va boshqa texnik ketma-ketliklarni topib olib tashlaydi.

```
ILLUMINACLIP:<fastaWithAdaptersEtc>:<seedMismatches>:<palindromeClipThreshold>:<simpleClipThreshold>[:<minAdapterLength>[:<keepBothReads>]]
```

#### Ikki bosqichli qidiruv (seed + full alignment)

1. **Seed:** adapterning qisqa bo'laklari (**maksimum 16 bp**) read'ning har bir mumkin bo'lgan o'rnida tekshiriladi. Bo'lak mukammal yoki yetarlicha yaqin (`seedMismatches` ruxsat doirasida) mos kelsa, keyingi bosqichga o'tiladi.
2. **Full alignment:** read va adapterning to'liq moslashuvi **ball** bilan baholanadi.

Bu ikki bosqichli yondashuv tezlikni oshiradi: seed hisobi juda tez, to'liq ball esa nisbatan kam hisoblanadi.

#### Ball formulasi

- Har bir **mos** asos: **+0.6**
- Har bir **mos kelmagan** asos: **−Q/10** (Q — shu asosning sifati)

Sifat hisobga olingani uchun sekvenator xatosidan kelib chiqqan nomuvofiqliklar kamroq jazolanadi.

| Moslashuv uzunligi (xatosiz) | Ball |
|---|---|
| 12 asos | ≈ 7.2 |
| 17 asos | ≈ 10.2 |
| 25 asos | 15 |
| ~50 asos | ≈ 30 |

Shuning uchun qo'llanma **simple** rejimda **7–15** chegarani tavsiya qiladi; **palindrome** rejimida (uzunroq moslashuv mumkin) chegara **~30** bo'lishi mumkin.

#### Parametrlar

| Parametr | Ma'nosi |
|---|---|
| `fastaWithAdaptersEtc` | Adapter va PCR ketma-ketliklari fayli. Fayldagi nomlar ketma-ketlik qanday ishlatilishini belgilaydi (6-bo'lim) |
| `seedMismatches` | Seed bosqichida ruxsat etilgan maksimal nomuvofiqlik soni |
| `palindromeClipThreshold` | PE palindrome moslashuvi qanchalik aniq bo'lishi kerakligi (ball chegarasi) |
| `simpleClipThreshold` | Har qanday adapter ketma-ketligining read bilan moslashuvi qanchalik aniq bo'lishi kerakligi (ball chegarasi) |
| `minAdapterLength` | **Faqat palindrome.** Topilgan adapterning minimal uzunligi. Standart: **8** (tarixiy sabab). Palindrome rejimida soxta musbat juda kam bo'lgani uchun 1 gacha kamaytirsa bo'ladi |
| `keepBothReads` | **Faqat palindrome.** Read-through topilib adapter olib tashlangach, reverse read forward read bilan bir xil ma'lumotni (reverse-complement) saqlaydi. Standart holatda reverse read **tashlanadi**. `true` bo'lsa, ikkalasi saqlanadi |

> **SE rejimida** `palindromeClipThreshold` hech qanday ta'sir qilmaydi, lekin sintaksis uchun baribir qiymat yozish kerak (misol: `2:30:10`).

#### Palindrome va Simple rejimlari

- **Simple:** har bir adapter ketma-ketligi read bilan tekshiriladi; yetarlicha aniq moslik topilsa, read tegishli joyidan kesiladi.
- **Palindrome:** «read-through» uchun mo'ljallangan. Adapterlar forward va reverse read'larning boshiga *in silico* «ulanadi», va ikkala ketma-ketlik (adapter + read) bir-biriga tekislanadi. Agar ular read-through ko'rsatuvchi tarzda mos kelsa, forward read kesiladi, reverse read esa tashlanadi (unda yangi ma'lumot yo'q).

Palindrome yondashuvi ancha ishonchli, chunki tekshiriladigan hudud uzun: ikkala adapter bir vaqtda va ikkala read'ning parcha qismlari ham tekshiriladi. Shu sababli adapterning **atigi bitta asosi** o'qilgan bo'lsa ham read-through aniqlanishi mumkin.

### 5.2. SLIDINGWINDOW

```
SLIDINGWINDOW:<windowSize>:<requiredQuality>
```

5' uchidan boshlab skanerlaydi va oyna ichidagi **o'rtacha** sifat chegaradan tushganda read'ni kesadi. Bir nechta asos o'rtachalanishi tufayli, bitta yomon asos keyingi sifatli ma'lumotni olib tashlashga olib kelmaydi.

Masalan, `SLIDINGWINDOW:4:15` — 4 asosli oyna, o'rtacha sifat 15 dan past tushsa kes.

### 5.3. MAXINFO

```
MAXINFO:<targetLength>:<strictness>
```

Moslashuvchan (adaptive) sifat bo'yicha kesish: uzun read saqlash foydasi bilan xatoli asoslarni saqlash narxini muvozanatlaydi. Read «qiymati» uchta omilga bog'liq:

1. **Minimal uzunlik:** read ketma-ketlikda noyob joylashishi uchun yetarlicha uzun bo'lishi kerak (odatda ~40 asos).
2. **Qo'shimcha uzunlik:** noyob joylashish uchun yetarli bo'lgandan keyin ham qo'shimcha asoslar foydali (assembly va variant qidirishda).
3. **Xatoga sezgirlik:** ba'zi vositalar bitta xatoga ham sezgir, boshqalari ko'p xatoga chidamli.

Parametrlar:

| Parametr | Ma'nosi |
|---|---|
| `targetLength` | Read joylashishini aniqlash uchun yetarli uzunlik |
| `strictness` | 0 dan 1 gacha. **Past (<0.2)** — uzunroq read'lar; **yuqori (>0.8)** — to'g'ri asoslar afzal |

Kesish read'ning 3' uchida, har bir mumkin pozitsiyada uch omilni birlashtirib hisoblanadi va eng yaxshi ball tanlanadi.

### 5.4. LEADING

```
LEADING:<quality>
```

Read **boshidan** sifati chegaradan past asoslarni olib tashlaydi: asos chegaradan past bo'lsa olib tashlanadi, keyingisi tekshiriladi, to'xtaguncha.

### 5.5. TRAILING

```
TRAILING:<quality>
```

Read **oxiridan** shunday qiladi. Qo'llanma: bu bilan Illumina'ning maxsus «low quality segment» (Q=2) hududlarini olib tashlash mumkin, **lekin SLIDINGWINDOW yoki MAXINFO tavsiya etiladi**.

### 5.6. CROP

```
CROP:<length>
```

Sifatdan qat'i nazar, read oxiridan asoslarni kesib, uzunligini **ko'pi bilan** `length` ga keltiradi (boshidan `length` ta asos saqlanadi). Undan keyingi bosqichlar read'ni yana qisqartirishi mumkin.

### 5.7. HEADCROP

```
HEADCROP:<length>
```

Sifatdan qat'i nazar, read **boshidan** `length` ta asosni olib tashlaydi. (Masalan, *Per base sequence content* da boshidagi notekislik uchun ishlatish mumkin.)

### 5.8. MINLEN

```
MINLEN:<length>
```

`length` dan qisqa read'larni tashlaydi. Kerak bo'lsa, **barcha boshqa bosqichlardan keyin** turishi kerak. Tashlangan read'lar Trimmomatic xulosasidagi «dropped reads» sonida hisoblanadi.

### 5.9. AVGQUAL

```
AVGQUAL:<quality>
```

O'rtacha sifati belgilangan darajadan past read'ni tashlaydi. (0.30 versiyada qo'shilgan.)

### 5.10. TOPHRED33 / TOPHRED64

Sifat kodlashini phred+33 yoki phred+64 ga o'zgartiradi. Qo'shimcha parametr yo'q.

### 5.11. Bosqichlar jadvali

| Bosqich | Qayerdan | Nima qiladi | Read'ni tashlaydimi? |
|---|---|---|---|
| ILLUMINACLIP | 3' (asosan) | Adapterni kesadi | Yo'q (PE palindromda reverse read tashlanishi mumkin) |
| SLIDINGWINDOW | 5' dan skanerlab, 3' ni kesadi | O'rtacha sifat bo'yicha | Yo'q |
| MAXINFO | 3' | Uzunlik va xato muvozanati | Yo'q |
| LEADING | 5' | Past sifatli asoslar | Yo'q |
| TRAILING | 3' | Past sifatli asoslar | Yo'q |
| CROP | 3' | Uzunlikka qirqadi | Yo'q |
| HEADCROP | 5' | N ta asosni olib tashlaydi | Yo'q |
| MINLEN | — | Qisqa read'ni tashlaydi | **Ha** |
| AVGQUAL | — | O'rtacha sifat past read'ni tashlaydi | **Ha** |
| TOPHRED33/64 | — | Kodlash o'zgartiradi | Yo'q |

---

## 6. Adapter fayllari

### 6.1. Tayyor fayllar

Illumina adapterlari mualliflik huquqi bilan himoyalangan, lekin Trimmomatic bilan tarqatishga ruxsat olingan.

| Fayl | Qaysi kutubxona |
|---|---|
| `TruSeq2-SE.fa`, `TruSeq2-PE.fa` | TruSeq2 (GAII mashinalari) |
| `TruSeq3-SE.fa`, `TruSeq3-PE.fa` | TruSeq3 (HiSeq va MiSeq) |

Qo'llanma qoidasi: **yangiroq kutubxonalar TruSeq3 ishlatadi**, lekin bu xizmat ko'rsatuvchiga bog'liq.

### 6.2. FastQC bilan to'g'ri faylni tanlash

FastQC'ning *Overrepresented sequences* hisoboti qaysi fayl mosligini ko'rsatadi:

| FastQC «Possible Source» | Kerakli fayl |
|---|---|
| «Illumina Single End» / «Illumina Paired End» | `TruSeq2-SE.fa` / `TruSeq2-PE.fa` |
| «**TruSeq Universal Adapter**» yoki «**TruSeq Adapter, Index …**» | `TruSeq3-SE.fa` (single-end) yoki `TruSeq3-PE.fa` (paired-end) |
| «Illumina Multiplexing …» va RNA kutubxonalari | Qo'llanma versiyasida **kiritilmagan** |

> Qo'llanma Nextera uchun ham hozircha fayl yo'qligini aytadi. Yangiroq Trimmomatic versiyalarida `NexteraPE-PE.fa` ham keladi.

### 6.3. Maxsus adapter fayli yaratish

FASTA fayldagi **ketma-ketlik nomlari** uning ishlatilishini belgilaydi.

**Palindrome uchun** moslangan juft kerak:

- Nomlar `Prefix` bilan boshlanadi.
- Forward adapter `/1` bilan, reverse adapter `/2` bilan tugaydi.
- `Prefix` va `/1` (yoki `/2`) orasidagi qism juftlikda **aynan bir xil** bo'lishi kerak.

```
>PrefixPE/1
AGATGTGTATAAGAGACAG
>PrefixPE/2
AGATGTGTATAAGAGACAG
```

Juft hosil qilmaydigan `Prefix…` ketma-ketliklar hozircha simple rejimda ishlatiladi (kelajakda xato bo'ladi).

**Simple uchun:** palindromda ishlatilmagan barcha ketma-ketliklar.

- Nomi `/1` yoki `/2` bilan tugasa, faqat forward yoki faqat reverse read bilan tekshiriladi.
- Boshqa barcha ketma-ketliklar ikkala read bilan tekshiriladi.

**Reverse-complement avtomatik tekshirilmaydi.** Kerak bo'lsa, uni alohida nom bilan faylga qo'shish kerak. Qo'llanmadagi `TruSeq2-PE.fa` misoli:

```
>PCR_Primer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR_Primer1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
```

Sabab: ba'zi ketma-ketliklar bir yo'nalishda ikkinchisidan ancha ko'proq uchraydi.

---

## 7. Qo'llanmadagi misollar

### 7.1. Paired-end

```bash
java -jar trimmomatic-0.30.jar PE s_1_1_sequence.txt.gz s_1_2_sequence.txt.gz \
    lane1_forward_paired.fq.gz lane1_forward_unpaired.fq.gz \
    lane1_reverse_paired.fq.gz lane1_reverse_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Bajarilish tartibi:

1. **ILLUMINACLIP** (`TruSeq3-PE.fa`): avval 16 asosli seed, ko'pi bilan **2** nomuvofiqlik bilan. Seed uzaytiriladi va PE read'lar uchun ball **30** (~50 asos), single-end uchun **10** (~17 asos) ga yetsa, adapter kesiladi.
2. **LEADING:3** — boshidagi sifati 3 dan past (yoki N) asoslarni olib tashlaydi.
3. **TRAILING:3** — oxiridagi sifati 3 dan past (yoki N) asoslarni olib tashlaydi.
4. **SLIDINGWINDOW:4:15** — 4 asosli oyna, o'rtacha sifat 15 dan tushsa kesadi.
5. **MINLEN:36** — shu bosqichlardan keyin 36 asosdan qisqa read'larni tashlaydi.

### 7.2. Single-end

```bash
java -jar trimmomatic-0.30.jar SE s_1_1_sequence.txt.gz lane1_forward.fq.gz \
    ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Xuddi shu bosqichlar, lekin single-end adapter fayli bilan. `:30:` parametrining ta'siri yo'q, lekin qiymat berilishi shart.

---

## 8. Algoritmlar: pseudocode, murakkablik, kod

> **Muhim eslatma.** Quyidagi kodlar — Trimmomatic g'oyasini tushunish uchun **o'quv modeli**. Ular dasturning to'liq nusxasi emas: haqiqiy Trimmomatic'da ILLUMINACLIP'dagi seed bosqichi, PE palindrome, SLIDINGWINDOW'ning read oxiridagi chegara holatlari (0.30 versiya o'zgarishlarida shunga tegishli tuzatish bor) va boshqa detallar bor. Natijalar bir necha asosga farq qilishi mumkin. Loyiha natijalari uchun doim haqiqiy Trimmomatic'dan foydalan.

### 8.1. LEADING va TRAILING

**Pseudocode:**

```
FUNCTION Leading(seq, q, threshold):
    i ← 0
    WHILE i < LENGTH(q) AND q[i] < threshold:
        i ← i + 1
    RETURN seq[i..], q[i..]

FUNCTION Trailing(seq, q, threshold):
    j ← LENGTH(q)
    WHILE j > 0 AND q[j-1] < threshold:
        j ← j - 1
    RETURN seq[..j), q[..j)
```

**Murakkablik:** vaqt **O(k)**, bunda k — kesilgan asoslar soni (eng yomon holda O(n)). Qo'shimcha xotira O(1) (indekslar bilan ishlansa).

### 8.2. SLIDINGWINDOW

**G'oya:** oynaning yig'indisini har safar qayta hisoblash o'rniga, oyna siljiganda bitta asosni qo'shib, bittasini ayiramiz (rolling sum).

```
FUNCTION SlidingWindow(seq, q, w, threshold):
    n ← LENGTH(q)
    IF n < w:
        RETURN (butun read o'rtachasi ≥ threshold) ? read : bo'sh
    window ← SUM(q[0 .. w-1])
    FOR i ← 0 TO n - w:
        IF i > 0:
            window ← window + q[i + w - 1] - q[i - 1]
        IF window < threshold × w:        // o'rtacha < threshold
            RETURN seq[0 .. i)            // oyna boshigacha saqla
    RETURN seq
```

O'rtachani butun sonlarda saqlash uchun `average < thr` ni `sum < thr × w` ko'rinishida tekshiramiz (bo'lishsiz, aniqroq).

**Murakkablik:**

- Naiv usul (har oynada qayta yig'ish): **O(n · w)**.
- Rolling sum: **O(n)** vaqt, **O(1)** qo'shimcha xotira.
- Butun fayl uchun N ta read: **O(N · n)**. Sening faylingda N = 2·10⁶, n = 150 → ~3·10⁸ elementar amal.

**Qo'lda hisoblash misoli.** `SLIDINGWINDOW:4:15`, sifatlar: `40 40 40 38 12 10 8 5 3 2`

| i | Oyna | Yig'indi | O'rtacha | Natija |
|---|---|---|---|---|
| 0 | 40 40 40 38 | 158 | 39.5 | davom |
| 1 | 40 40 38 12 | 130 | 32.5 | davom |
| 2 | 40 38 12 10 | 100 | 25.0 | davom |
| 3 | 38 12 10 8 | 68 | 17.0 | davom |
| 4 | 12 10 8 5 | 35 | 8.75 | **< 15 → kes** |

Natija: dastlabki **4** asos saqlanadi (`i = 4` gacha).

### 8.3. MINLEN

```
FUNCTION MinLen(seq, L):
    IF LENGTH(seq) < L: RETURN DROP
    RETURN seq
```

Murakkablik O(1) (uzunlik tayyor saqlangan bo'lsa).

### 8.4. ILLUMINACLIP (simple rejim, soddalashtirilgan)

```
FUNCTION ClipAdapter(read, q, adapter, threshold):
    n ← LENGTH(read);  m ← LENGTH(adapter)
    bestScore ← -∞;  bestPos ← NONE
    FOR i ← 0 TO n - 1:
        L ← MIN(m, n - i)                  // 3' uchida adapter qisman sig'adi
        score ← 0
        FOR k ← 0 TO L - 1:
            IF read[i + k] = adapter[k]: score ← score + 0.6
            ELSE:                        score ← score - q[i + k] / 10
        IF score ≥ threshold AND score > bestScore:
            bestScore ← score;  bestPos ← i
    IF bestPos = NONE: RETURN read
    RETURN read[0 .. bestPos)
```

**Murakkablik:**

- O'quv modeli: **O(n · m)** (n — read, m — adapter uzunligi).
- Haqiqiy Trimmomatic: avval seed (≤16 bp) tekshiruvi har pozitsiyada O(16), to'liq ball esa faqat seed o'tgan pozitsiyalarda. Shuning uchun amalda tezroq.

**Demo natijasi (threshold = 10, hamma asoslar Q=40):**

| Read oxiridagi adapter qismi | Ball | Natija |
|---|---|---|
| 25 asos (xatosiz) | 25 × 0.6 = **15.0** | kesiladi (pozitsiya 40) |
| 12 asos (xatosiz) | 12 × 0.6 = **7.2** | **kesilmaydi** (chegaradan past) |

Bu — 5.1-dagi sezgirlik va o'ziga xoslik kelishuvining amaldagi ko'rinishi: qisqa adapter bo'laklari (12 asos) 10-chegarada qolib ketadi. Chegarani 7 ga tushirsang, ular topiladi, lekin soxta moslashuvlar ham ko'payadi.

### 8.5. Kodlar: LEADING + TRAILING + SLIDINGWINDOW + MINLEN

Kirish: `sequence<TAB>quality` (phred+33) qatorlari; chiqish: kesilgan ketma-ketlik yoki `DROPPED`.

#### Python

```python
#!/usr/bin/env python3
"""Trimmomatic'ning asosiy bosqichlarining o'quv (soddalashtirilgan) modeli.
Kirish: TSV (sequence <TAB> quality[Phred+33]) -> chiqish: kesilgan sequence yoki DROPPED.
"""
import sys

def phred(qstr, offset=33):
    return [ord(c) - offset for c in qstr]

def leading(seq, q, thr):
    i = 0
    while i < len(q) and q[i] < thr:
        i += 1
    return seq[i:], q[i:]

def trailing(seq, q, thr):
    j = len(q)
    while j > 0 and q[j - 1] < thr:
        j -= 1
    return seq[:j], q[:j]

def sliding_window(seq, q, w, thr):
    n = len(q)
    if n < w:                      # read oynadan qisqa: butun read bitta oyna
        return (seq, q) if n and sum(q) / n >= thr else ("", [])
    window = sum(q[:w])
    for i in range(n - w + 1):
        if i > 0:
            window += q[i + w - 1] - q[i - 1]
        if window < thr * w:       # o'rtacha < thr  <=>  yig'indi < thr*w
            return seq[:i], q[:i]
    return seq, q

def minlen(seq, q, length):
    return (seq, q) if len(seq) >= length else None

def run(seq, qstr, steps):
    q = phred(qstr)
    for name, *args in steps:
        if name == "LEADING":
            seq, q = leading(seq, q, int(args[0]))
        elif name == "TRAILING":
            seq, q = trailing(seq, q, int(args[0]))
        elif name == "SLIDINGWINDOW":
            seq, q = sliding_window(seq, q, int(args[0]), int(args[1]))
        elif name == "MINLEN":
            if minlen(seq, q, int(args[0])) is None:
                return None
    return seq

if __name__ == "__main__":
    steps = [s.split(":") for s in sys.argv[1:]]
    for line in sys.stdin:
        s, qs = line.rstrip("\n").split("\t")
        r = run(s, qs, steps)
        print("DROPPED" if r is None else r)
```

#### C++

```cpp
// Trimmomatic bosqichlarining o'quv modeli (C++17). Kirish: stdin'da "sequence<TAB>quality".
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Read { std::string seq; std::vector<int> q; };

static void leading(Read& r, int thr) {
    size_t i = 0;
    while (i < r.q.size() && r.q[i] < thr) ++i;
    r.seq.erase(0, i); r.q.erase(r.q.begin(), r.q.begin() + i);
}
static void trailing(Read& r, int thr) {
    size_t j = r.q.size();
    while (j > 0 && r.q[j - 1] < thr) --j;
    r.seq.resize(j); r.q.resize(j);
}
static void sliding_window(Read& r, int w, int thr) {
    int n = (int)r.q.size();
    if (n < w) {
        long s = 0; for (int x : r.q) s += x;
        if (n == 0 || (double)s / n < thr) { r.seq.clear(); r.q.clear(); }
        return;
    }
    long window = 0;
    for (int i = 0; i < w; ++i) window += r.q[i];
    for (int i = 0; i + w <= n; ++i) {
        if (i > 0) window += r.q[i + w - 1] - r.q[i - 1];
        if (window < (long)thr * w) { r.seq.resize(i); r.q.resize(i); return; }
    }
}

int main(int argc, char* argv[]) {
    std::vector<std::vector<std::string>> steps;
    for (int i = 1; i < argc; ++i) {
        std::vector<std::string> parts; std::stringstream ss(argv[i]); std::string p;
        while (std::getline(ss, p, ':')) parts.push_back(p);
        steps.push_back(parts);
    }
    std::string line;
    while (std::getline(std::cin, line)) {
        size_t tab = line.find('\t');
        Read r; r.seq = line.substr(0, tab);
        for (char c : line.substr(tab + 1)) r.q.push_back(c - 33);
        bool dropped = false;
        for (auto& s : steps) {
            if (s[0] == "LEADING") leading(r, std::stoi(s[1]));
            else if (s[0] == "TRAILING") trailing(r, std::stoi(s[1]));
            else if (s[0] == "SLIDINGWINDOW") sliding_window(r, std::stoi(s[1]), std::stoi(s[2]));
            else if (s[0] == "MINLEN") { if ((int)r.seq.size() < std::stoi(s[1])) { dropped = true; break; } }
        }
        std::cout << (dropped ? "DROPPED" : r.seq) << "\n";
    }
}
```

#### C

```c
/* Trimmomatic bosqichlarining o'quv modeli (C). Kirish: stdin'da "sequence<TAB>quality". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXL 100000

int main(int argc, char *argv[]) {
    static char line[2 * MAXL + 2];
    static int q[MAXL];
    while (fgets(line, sizeof line, stdin)) {
        line[strcspn(line, "\n")] = '\0';
        char *tab = strchr(line, '\t');
        if (!tab) continue;
        *tab = '\0';
        char *seq = line, *qs = tab + 1;
        int n = (int)strlen(seq);
        for (int i = 0; i < n; ++i) q[i] = qs[i] - 33;
        int start = 0, end = n, dropped = 0;           /* saqlanadigan qism: [start, end) */
        for (int a = 1; a < argc && !dropped; ++a) {
            char step[64]; strncpy(step, argv[a], 63); step[63] = '\0';
            char *name = strtok(step, ":");
            char *p1 = strtok(NULL, ":");
            char *p2 = strtok(NULL, ":");
            if (!strcmp(name, "LEADING")) {
                int thr = atoi(p1);
                while (start < end && q[start] < thr) ++start;
            } else if (!strcmp(name, "TRAILING")) {
                int thr = atoi(p1);
                while (end > start && q[end - 1] < thr) --end;
            } else if (!strcmp(name, "SLIDINGWINDOW")) {
                int w = atoi(p1), thr = atoi(p2), len = end - start;
                if (len < w) {
                    long s = 0; for (int i = start; i < end; ++i) s += q[i];
                    if (len == 0 || (double)s / len < thr) end = start;
                } else {
                    long window = 0;
                    for (int i = 0; i < w; ++i) window += q[start + i];
                    for (int i = 0; i + w <= len; ++i) {
                        if (i > 0) window += q[start + i + w - 1] - q[start + i - 1];
                        if (window < (long)thr * w) { end = start + i; break; }
                    }
                }
            } else if (!strcmp(name, "MINLEN")) {
                if (end - start < atoi(p1)) dropped = 1;
            }
        }
        if (dropped) puts("DROPPED");
        else printf("%.*s\n", end - start, seq + start);
    }
    return 0;
}
```

**Sinov** (uchala dastur bir xil natija berdi):

```bash
python3 trim.py LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:10 < reads.tsv
g++ -std=c++17 -O2 -o trim_cpp trim.cpp && ./trim_cpp LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:10 < reads.tsv
gcc -O2 -o trim_c trim.c && ./trim_c LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:10 < reads.tsv
```

Test ma'lumotlarning tavsifi va natijasi:

| # | Read | Sifat tuzilishi | Natija | Sabab |
|---|---|---|---|---|
| 1 | 50 nt | boshida `##` (Q=2), oxirida 10 ta `#` | **38 nt** | LEADING 2 ta, TRAILING 10 ta asosni oldi |
| 2 | 30 nt | hammasi `I` (Q=40) | **30 nt** | o'zgarmadi |
| 3 | 20 nt | 8 ta Q=40, keyin 12 ta Q=10 | **DROPPED** | SLIDINGWINDOW 8-asosdan kesdi → 8 < MINLEN:10 |
| 4 | 12 nt | hammasi `#` (Q=2) | **DROPPED** | LEADING hammasini olib tashladi |

### 8.6. Kodlar: ILLUMINACLIP (simple, soddalashtirilgan)

#### Python

```python
#!/usr/bin/env python3
"""ILLUMINACLIP 'simple' rejimining o'quv modeli (seed bosqichisiz).
Ball: mos kelgan asos +0.6, mos kelmagan asos -Q/10.
"""
def clip_adapter(read, qual, adapter, threshold=10.0):
    n, m = len(read), len(adapter)
    best_score, best_pos = float("-inf"), None
    for i in range(n):                       # adapter boshlanishi uchun nomzod pozitsiya
        L = min(m, n - i)                    # 3' uchida adapter qisman sig'adi
        score = 0.0
        for k in range(L):
            score += 0.6 if read[i + k] == adapter[k] else -qual[i + k] / 10
        if score >= threshold and score > best_score:
            best_score, best_pos = score, i
    return (read[:best_pos], best_pos, best_score) if best_pos is not None else (read, None, best_score)

if __name__ == "__main__":
    import random
    random.seed(7)
    ad = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    frag = "".join(random.choice("ACGT") for _ in range(40))
    for tail in (25, 12):
        read = frag + ad[:tail]
        out, pos, sc = clip_adapter(read, [40] * len(read), ad, 10.0)
        print(f"adapter qismi={tail} nt, read={len(read)} nt -> pos={pos}, ball={sc:.1f}, natija uzunligi={len(out)}")
```

#### C++

```cpp
// ILLUMINACLIP 'simple' rejimining o'quv modeli (C++17), seed bosqichisiz.
#include <iostream>
#include <limits>
#include <string>
#include <vector>

int clip_pos(const std::string& read, const std::vector<int>& q, const std::string& ad, double thr) {
    int n = (int)read.size(), m = (int)ad.size(), best = -1;
    double best_score = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
        int L = std::min(m, n - i);
        double s = 0;
        for (int k = 0; k < L; ++k) s += (read[i + k] == ad[k]) ? 0.6 : -q[i + k] / 10.0;
        if (s >= thr && s > best_score) { best_score = s; best = i; }
    }
    return best;   // -1: adapter topilmadi
}

int main() {
    std::string ad = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
    std::string frag = "ACGTTGCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGA";   // 40 nt
    for (int tail : {25, 12}) {
        std::string read = frag + ad.substr(0, tail);
        std::vector<int> q(read.size(), 40);
        int pos = clip_pos(read, q, ad, 10.0);
        std::cout << "adapter qismi=" << tail << " nt -> pos=" << pos << "\n";
    }
}
```

---

## 9. Sening loyihang uchun amaliy qism

FastQC hisobotingdagi faktlar va ularga mos Trimmomatic qarorlari:

| FastQC topilmasi | Qo'llanmadagi mos vosita |
|---|---|
| Adapter Content FAIL (Universal Adapter ~83%); Overrepresented'da «**TruSeq Adapter, Index 2**» | Qo'llanma jadvaliga ko'ra bu **TruSeq3** kutubxonasi → `TruSeq3-SE.fa` (yoki PE bo'lsa `TruSeq3-PE.fa`) |
| Oxirgi pozitsiyalarda Q=2 (10-persentil) | Illumina «low quality segment». **TRAILING:3** olib tashlaydi, lekin qo'llanma **SLIDINGWINDOW** ni afzal ko'radi — ikkalasini birga ishlatish xavfsiz |
| 110-pozitsiyadan keyin sifat tushishi | **SLIDINGWINDOW:4:15** |
| Boshidagi (1–10) notekis tarkib | Kesish bilan to'liq tuzalmaydi; ixtiyoriy **HEADCROP** bilan sinab ko'rish mumkin (lekin o'zgartirish sababini P2P'da tushuntira olishing kerak) |

### 9.1. Single-end buyrug'i

```bash
# 2.2: faqat ILLUMINACLIP
trimmomatic SE data/SRR3341836_1.fastq results/2.2_trimmed.fastq \
    ILLUMINACLIP:data/TruSeq3-SE.fa:2:30:10

# 2.7: to'liq (qo'llanma namunasi)
trimmomatic SE data/SRR3341836_1.fastq results/2.7_trimmed.fastq \
    ILLUMINACLIP:data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Diskni tejash uchun (4.1-bo'lim) chiqish nomlarini `.fastq.gz` qilsang bo'ladi:

```bash
trimmomatic SE data/SRR3341836_1.fastq results/2.7_trimmed.fastq.gz \
    ILLUMINACLIP:data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 9.2. Agar ma'lumot paired-end bo'lsa

Fayl nomi `_1` bilan tugagani PE bo'lishi mumkinligini bildiradi (NCBI sahifasidagi **Layout** maydonini tekshir). Unda:

```bash
trimmomatic PE -baseout results/2.7.fastq.gz \
    data/SRR3341836_1.fastq data/SRR3341836_2.fastq \
    ILLUMINACLIP:data/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Natijada `2.7_1P`, `2.7_1U`, `2.7_2P`, `2.7_2U` fayllari hosil bo'ladi.

### 9.3. Natijani qayerdan o'qish

Trimmomatic oxirida xulosa chiqaradi (`2> results/2.2.txt` bilan faylga saqlanadi):

```
Input Reads: ...  Surviving: ... (..%)  Dropped: ... (..%)
```

P2P'da shu raqamlar «hisobot yaxshilandi, lekin X% read yo'qotildi» degan javobning dalili bo'ladi.

---

## 10. Mashqlar va javoblar

**1-mashq.** Sifatlar `40 40 40 38 12 10 8 5 3 2`, bosqich `SLIDINGWINDOW:4:15`. Nechta asos saqlanadi?

**2-mashq.** `5`, `I`, `#`, `+` belgilarining Phred+33 sifatini top.

**3-mashq.** ILLUMINACLIP simple rejimida chegara 10. Xatosiz moslashuvda kamida nechta asos kerak? 15 asos yetadimi? 20 asos-chi?

**4-mashq.** Adapter bilan moslashuvda 1 ta nomuvofiqlik bor (Q=40 asosda). Chegara 10 bo'lsa, kamida nechta mos asos kerak?

**5-mashq.** Nima uchun `MINLEN` oxirida turishi kerak? ILLUMINACLIP nima uchun boshida?

**6-mashq.** Nima uchun SE rejimida `ILLUMINACLIP:file:2:30:10` dagi `30` yozilishi kerak, garchi ishlatilmasa ham?

<details>
<summary>Javoblar</summary>

1. **4** asos (i=4 oynasi o'rtachasi 8.75 < 15; oldingi oynalar ≥ 17).
2. `5` → 20; `I` → 40; `#` → 2; `+` → 10.
3. Ball ≥ 10 uchun 0.6·m ≥ 10 → m ≥ 16.67 → **kamida 17 asos** (17×0.6 = 10.2). 15 asos (9.0) **yetmaydi**; 20 asos (12.0) **yetadi**.
4. 0.6·m − 4 ≥ 10 → m ≥ 23.33 → **24 ta mos asos** (24×0.6 − 4 = 10.4; 23 ta bo'lsa 9.8 — yetmaydi).
5. MINLEN uzunlikni *boshqa bosqichlar read'ni qisqartirgandan keyin* tekshiradi; boshida tursa, keyin qisqarib qolgan read'lar tekshiruvdan o'tib ketadi. ILLUMINACLIP boshida, chunki sifat bo'yicha kesish adapterning bir qismini olib tashlasa, qolgan qisman adapterni topish qiyinlashadi.
6. Chunki ILLUMINACLIP sintaksisi uch ketma-ket parametrni (`seedMismatches:palindromeClipThreshold:simpleClipThreshold`) talab qiladi; palindrome chegarasi faqat PE'da ishlatiladi, SE'da ta'siri yo'q, lekin pozitsiya uchun qiymat kerak.

</details>

---

## 11. Qisqa ma'lumotnoma

```
SE:  trimmomatic SE [-phred33] [-trimlog f] in out STEP...
PE:  trimmomatic PE [-phred33] in1 in2 out1P out1U out2P out2U STEP...
     (yoki -baseout nom.fq.gz)

ILLUMINACLIP:fasta:seedMM:palindromeThr:simpleThr[:minAdapterLen[:keepBoth]]
SLIDINGWINDOW:oyna:sifat        MAXINFO:targetLen:strictness (0..1)
LEADING:sifat   TRAILING:sifat  CROP:uzunlik   HEADCROP:uzunlik
MINLEN:uzunlik  AVGQUAL:sifat   TOPHRED33   TOPHRED64

Standart (qo'llanma): ILLUMINACLIP:TruSeq3-xx.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
Ball: mos +0.6, nomuvofiqlik -Q/10 | simple chegara 7..15 | palindrome ~30
Tartib: ILLUMINACLIP birinchi, MINLEN oxirgi
```

---

## 12. Versiyalar tarixi

| Versiya | O'zgarishlar |
|---|---|
| 0.32 | Sifat kodlashi avtomatik aniqlanadi (avval standart `phred64` edi) |
| 0.30 | AVGQUAL qo'shildi; SLIDINGWINDOW'da read oxiridagi «yarim oyna» kesish xatosi tuzatildi (0.27 da kiritilgan); ILLUMINACLIP simple moslashuvi tezlashtirildi |
| 0.27 | Palindrome'ga «keep reverse» va «minimum adapter length» qo'shildi; `-jar` bilan ishga tushirish; MAXINFO; simple rejimda qisqa adapterlar yaxshilandi; FASTQ parser tezlashtirildi |
| 0.25 | Birlashtirilgan (concatenated) gzip, kirish fayllarini faqat bir marta o'qish, ko'p oqimli rejimda xato xabari |
| 0.22 | ZIP va BZIP2 qo'llab-quvvatlash; simple rejimda qisqa adapterlar uchun asosiy qo'llab-quvvatlash |
| 0.20 | Ko'p oqimlilik, qisqacha kesish statistikasi, ILLUMINACLIP ball xatosi tuzatildi |
