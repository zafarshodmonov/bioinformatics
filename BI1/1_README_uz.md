# pUC19 plazmidasini tahlil qilish

Bu — dasturdagi birinchi loyihang. Ish davomida sen genomik bioinformatik nima bilan shug'ullanishi bilan tanishasan, sekvenirlash (sequencing) natijasida olinadigan ma'lumotlarning asosiy formatlari hamda sifat nazorati (quality control) va filtrlash tamoyillarini o'zlashtirasan. Sifatni tahlil qilish uchun FastQC, filtrlash uchun Trimmomatic, ketma-ketliklarni qayta ishlash (GC-tarkibni aniqlash) uchun esa Biopython'dan foydalanishni o'rganasan.

💡 [Bu yerni bos](https://new.oprosso.net/p/4cb31ec3f47a4596bc758ea1861fb624), **shu loyiha bo'yicha bizga fikr-mulohaza bildirish uchun**. So'rovnoma anonim va jamoamizga o'qitishni yaxshilashda yordam beradi. Loyihani tugatgan zahoti to'ldirishni tavsiya qilamiz.

## Mundarija

  - [«21-Maktab»da qanday o'qish kerak](#21-maktabda-qanday-oqish-kerak)
  - [Loyiha bilan qanday ishlash kerak](#loyiha-bilan-qanday-ishlash-kerak)
  - [Chapter I](#chapter-i)
    - [Bioinformatikning bir kuni](#bioinformatikning-bir-kuni)
  - [Chapter II](#chapter-ii)
    - [**1-topshiriq.** Sekvenirlash sifatini baholash](#1-topshiriq-sekvenirlash-sifatini-baholash)
    - [**2-topshiriq.** Sekvenirlash sifatini yaxshilash (reads'larni tozalash)](#2-topshiriq-sekvenirlash-sifatini-yaxshilash-readslarni-tozalash)
    - [**3-topshiriq.** GC-tarkibni baholash](#3-topshiriq-gc-tarkibni-baholash)


## «21-Maktab»da qanday o'qish kerak

- Bu yerda seni katta erkinlikka ega noyob ta'lim tajribasi kutmoqda. Sen topshiriq olasan va yechim yo'llarini o'zing izlaysan: qulay bo'lgan har qanday usuldan — internet resurslaridan yoki neyrotarmoqlardan (masalan, GigaChat) foydalanishing mumkin. Biroq ma'lumot sifatiga e'tiborli bo'l: tekshir, o'yla, tahlil qil, solishtir.
- Peer-to-Peer (P2P) o'zaro ta'lim — boshqa peer'lar bilan bilim va tajriba almashinuvi bo'lib, bunda har kim ham o'qituvchi, ham o'quvchi bo'ladi. Bunday yondashuv bir-biringdan o'rganib, materialni chuqurroq tushunishga imkon beradi.
- O'zingni erkin his qil va yordam so'rashdan tortinma — atrofingdagilar ham bu yo'ldan birinchi marta o'tmoqda. Tajriba va g'oyalaringni boshqalar bilan ulash. Jamiyatimizdagi barcha yangiliklardan xabardor bo'lish uchun RocketChat'ga qo'shil.
- Agar birovning yechimini ko'chirsang, o'qishingning hech qanday ma'nosi qolmaydi. Boshqalarning yordamidan foydalansang ham, doim nima uchun, qanday va nima maqsadda ekanini oxirigacha tushunib ol. Xato qilishdan qo'rqma.
- Topshiriq bajarib bo'lmaydigandek tuyuladimi? Tanaffus qil, havo almashtir, boshingni «qayta yuklash» — bu ko'pchilikka yordam bergan. Ehtimol, yechim o'zi keladi.
- Faqat o'qish natijasi emas, jarayonning o'zi ham muhim. Masalani shunchaki yechish emas, uni QANDAY yechishni tushunish kerak.

## Loyiha bilan qanday ishlash kerak

- Bajarishdan oldin loyihani GitLab'dan xuddi shu nomdagi repozitoriyga clone qilish kerak.
- Barcha fayllar clone qilingan repozitoriyning _src/_ papkasida yaratilishi kerak.
- Loyihani clone qilgandan so'ng `develop` branch'ini yaratish va ishlab chiqishni shu branch'da olib borish kerak. Shundan keyin GitLab'ga ham `develop` branch'ini push qilish kerak.
- Direktoriyangda topshiriqlarda ko'rsatilganlardan boshqa fayllar bo'lmasligi kerak.

## Chapter I

### Bioinformatikning bir kuni

**9:45. Ertalab**

Budilnikning o'tkir ovozi meni uyqudan qayta-qayta uzib oladi — tushimda men A, T, G, C harflaridan iborat qochib borayotgan ketma-ketlikning orqasidan quvayotgan edim. Stolda kechagi sovib qolgan qahva, pochtada esa molekulyar biologlar jamoasidan o'qilmagan xat «yonib» turibdi: «SHOSHILINCH! pUC19 bilan eksperiment tugadi, oq koloniyalar tanlab olindi, plazmidalar ajratildi. Insert hosil bo'ldimi? Mutatsiyalar bormi? Yordam kerak! Ma'lumotlarni juda-juda kutamiz!».

Ular haftalab mehnat qilishdi, endi barcha umid — menda va mening algoritmlarimda. Men — ularning oxirgi instansiyasi, genomika davrining raqamli Sherlok Xolmsi.

Mening vazifam: xom, betartib ma'lumotlarni aniq biologik haqiqatga aylantirish. Vaqt ketmoqda. Kechagi qahva krujkasiga ko'zim tushadi. Aynan! Menga aynan shu kerak. Yangi qahva damlayman va nimadan boshlashni o'ylayman.

Birinchi navbatda — kunlik vazifalarni tizimlashtiraman. Asosiy algoritmlar va ishlash tamoyillarini esga tushiraman.

_P.S. Bioinformatik vazifalari haqida batafsil — Materials papkasida._

**10:00. Raqamli laboratoriya. «pUC19» ishini qabul qilish**

Tizimga kiraman. Ekranda — hozirgina kelgan gigabaytlab FASTQ fayllar.

Bu — «xom ma'lumotlar», sekvenator ishining raqamli aks-sadolari. Tasavvur qil: bakteriyalarning yashash bo'yicha kichik plazmida «ko'rsatmalar kitobi» (pUC19)ning minglab nusxalari millionlab mayda parchalar-gaplarga bo'lib tashlangan va skanerlangan.

Sekvenator (bizning «ma'lumot beruvchimiz») harakat qildi, lekin u mukammal emas. Bu parchalardagi har bir harf Phred Score (Q-score) bilan belgilangan — unga qanchalik ishonish mumkinligini aytuvchi raqam. Q=30? Yaxshi, xato ehtimoli atigi 1000 tadan 1 ta. Q=20? Yomonroq, 100 tadan 1 ta. Past Q-score — gumonli guvoh ko'rsatmalariga o'xshaydi: bu haqiqiy mutatsiyami yoki shunchaki «ma'lumot beruvchi»ning xatosimi?

Sifat nazorati — butun tergovning ishonchliligi poydevori. Yana meni ifloslanish (contamination) arvohi ta'qib qiladi — namunaga tasodifan tushib qolgan, soxta dalillarga o'xshash begona ketma-ketliklar. Umid qilamanki, barkodlar (har bir namunadagi noyob DNK-belgilar) o'z ishini qilgan.

Men sekvenator ishi sifatini baholashim, maqsadli insert borligini tasdiqlashim va kerakli oqsil sintezlanishini bashorat qilishim kerak.

_P.S. Molekulyar biologiya asoslari va sekvenirlash tamoyillari Materials papkasida batafsil ko'rib chiqilgan._

**10:15. «Hodisa joyi»ni ko'zdan kechirish: FastQC**

Birinchi qadam — falokat ko'lamini tushunish.

FastQC'ni ishga tushiraman — «dalillar»ni dastlabki ko'zdan kechirish uchun mening asbobim. U qizil va to'q sariq punktlar «chinqirib» turgan hisobot chizadi. Adapter Content hadsiz yuqori! Bu hodisa joyidan begona barmoq izlarini topishga o'xshaydi — sekvenirlash uchun DNK'ga tikilgan, ammo ma'lumotlarda qolib ketgan sun'iy ketma-ketliklar (adapterlar) parchalari. Ular hamma narsaga xalaqit beradi!

Per Base Sequence Quality read'lar oxiriga kelib sifat tushishini ko'rsatadi — «ma'lumot beruvchi» (bizning sekvenator) o'qishning oxiriga kelib charchagan.

Per Base Sequence Content gumonli darajada notekis ko'rinadi. Overrepresented sequences — ehtimol, ifloslanish yoki PCR artefaktlari izlaridir?

Hamkasb-biologlar chatda asabiy yozishadi: «Ma'lumotlarimiz qalay?».

Javob beraman: «Texnik axlat bor. Tozalayapman. Xabardor qilib turaman». Ularga javob kecha kerak edi, lekin shoshqaloqlik — aniqlikning dushmani. Bioinformatik metodik bo'lishi kerak.

_FastQC: sekvenirlash ma'lumotlari (FASTQ) sifatini vizualizatsiya qilish dasturi._

_Statuslar:_

_PASS (yashil)_

_WARN (to'q sariq — e'tibor talab qiladi)_

_FAIL (qizil — kritik)_

_Quyidagi parametrlarni tekshiradi:_

| **Per Base Sequence Quality** | **Per Sequence Quality Scores** | **Per Base Sequence Content** | **Adapter Content** | **Overrepresented sequences** | **Va boshqalar** |
| --- | --- | --- | --- | --- | --- |
| Read'lardagi har bir pozitsiya bo'yicha o'rtacha sifat (Q-score). | Read'larning o'rtacha sifat bo'yicha taqsimoti. | Har bir pozitsiyada A, T, G, C ning foiz tarkibi. | Adapter ketma-ketliklari bor read'lar foizi. | G'ayrioddiy tez-tez uchraydigan ketma-ketliklar (mumkin bo'lgan ifloslanish yoki artefakt). | Per seq GC content, Per base N content, Sequence length dist., Duplication levels va h.k. |

**10:50. Dalillarni «tozalash»: Trimmomatic**

Trimmomatic'ga o'tish vaqti keldi — mening raqamli qaychilarim. Mening vazifam — read'larning past sifatli qismlarini ehtiyotkorlik bilan kesib tashlash va adapterlarni olib tashlash. Bu texnik artefaktlar natijalarni jiddiy buzishi mumkin, shuning uchun puxta tozalashsiz bo'lmaydi. FastQC hisobotidan menga qaysi ketma-ketliklar xalaqit berayotgani allaqachon ma'lum. Trimmomatic'ning adapters papkasida aynan o'sha muammoli adapterlar bor faylni topaman. Demak, ILLUMINACLIP parametri kerak.

Parametrlarni ma'lumotlarga qarab sozlayman:

- 2 ta mos kelmaslikka (Mismatches) ruxsat beraman.
- O'xshashlik chegarasi (Palindrome / Simple clip threshold) — 30.
- Adapterning minimal uzunligi (MinAdapterLength) — 10.

Ishlov berishni ishga tushiraman va Trimmomatic ishlayotgan paytda yangi qahvani ho'plab turaman. Tozalash sifati juda muhim: «iflos» ma'lumotlarda assembly (yig'ish) noto'g'ri xulosalarga olib keladi.

Biologlar probirkadagi xafa mushukli mem yuborishadi. Ularning dardini tushunaman.

![cat](misc/images/cat.png)

*Manba: rasm Kandinsky 3.1 sun'iy intellekt algoritmi tomonidan yaratilgan.*

[_Trimmomatic: FASTQ fayllarni dastlabki qayta ishlash dasturi._](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

_Asosiy funksiyalar:_

**_ILLUMINACLIP_** _(adapterlarni olib tashlaydi. Adapter ketma-ketliklari fayli va parametrlar (mismatches, thresholds) talab qilinadi)._

**_SLIDINGWINDOW_** _(read'larni sifat bo'yicha sirpanuvchi oyna (sliding window) orqali kesadi)._

**_MINLEN_** _(juda qisqa read'larni tashlab yuboradi)._

**_TRAILING_** _(read'larning oxirini asoslar sifati bo'yicha kesadi)._

_Natija:_

- _tozalangan FASTQ fayllar;_
- _yaxshilangan ma'lumotlar sifati;_
- _ma'lumotlar alignment/assembly uchun tayyor._

**FASTQ fayllar → Trimmomatic → Toza ma'lumotlar**

**11:40. Qayta ko'zdan kechirish**

«Tozalash»dan keyingi FastQC. Fayllar tozalandi! FastQC'ni yana ishga tushiraman.

Adapter Content'ga hayajon bilan qarayman — yashil! Adapterlar yengildi. Per Base Sequence Quality tekislandi — past sifatli uchlar kesib tashlangan.

Garchi hali hammasi mukammal bo'lmasa-da: Per Base Sequence Content hamon to'q sariq zonada, read'lar uzunligi esa o'zgardi! Hmm… Ammo asosiysi — yaqqol axlatni olib tashladik. Ko'nglim yengillashdi.

Biologlarga yozaman: «Adapterlar olib tashlandi, sifat yaxshilandi. Assembly va insert tahliliga o'tyapman».

«🎉» degan javob smayli yoqimli, lekin oldinda hali ish ko'p…

**12:00. pUC19 plazmidasiga e'tibor**

Endi aynan nima sekvenirlanganini aniqlash kerak. Diqqat markazida — bizning «ashyoviy dalil»imiz, sun'iy pUC19 plazmidasi. Bu — DNK fragmentlarini klonlash uchun yaratilgan «biologlar uchun fleshka».

Davom etishdan oldin pUC19'ning asosiy xususiyatlarini esga tushiraman:
- plazmida tuzilishi;
- ishlash tamoyillari;
- tahlilning standart usullari.

_P.S. pUC19 plazmidasi tuzilishining to'liq tavsifi Materials papkasida._

**12:30**

Biologlar oq koloniyalar bo'yicha ma'lumot yuborishdi. Demak, MCS'ga nimadir kiritilgan. Biroq esda tutish muhim: ko'k-oq selektsiya (blue-white selection) — faqat dastlabki saralash bosqichi. To'g'ri fragment kiritilganmi? To'liq kiritilganmi? PCR, fragmentatsiya yoki klonlash paytida xatolar (mutatsiyalar) bo'lmaganmi?

Aynan shuning uchun pUC19 sekvenirlashni o'tkazish kerak: shu bakteriyalardan ajratilgan plazmida namunamizning aniq ketma-ketligini o'qib, kutilgan ketma-ketlik bilan solishtirish.

Bu — konstruktsiyani keyingi eksperimentlarda ishlatishdan oldingi yakuniy sifat nazorati.

**13:00. GC-tarkibni ekspress-tahlil qilish**

Insert aniq ketma-ketligining resurs talab qiluvchi assembly'si bajarilayotgan paytda, men tez bilvosita test qilaman — GC-tarkib tahlili. Bu ko'rsatkich butun plazmidada yoki qiziqtirgan fragmentda guanin (G) va sitozin (C) ning foiz miqdorini aks ettiradi.

- NCBI bazasidan pUC19 referens ketma-ketligini yuklab olaman.
- Biopython bilan Python'da ketma-ketlik bo'ylab yuguradigan va G hamda C foizini hisoblaydigan tezkor skript yozaman.
- Ishga tushiraman.
- Natijalarni etalon 51% bilan solishtiraman.

Natijalar:

- Og'ish <1% — ketma-ketlik kutilganga mos.
- Sezilarli farq (>5%) — qo'shimcha tekshiruv talab qiladi.

Joriy ma'lumotlar me'yor doirasida — namuna sifati foydasiga yana bir band.

_P.S. GC-tahlil ahamiyati va uni hisoblash usullari haqida batafsil — Materials papkasida._

**13:45. Ketma-ketlikni yig'ish va tahlil qilish**

Tozalangan sekvenirlash ma'lumotlari ishlov berishga tayyor. «Pazl»ni yig'ish vaqti keldi. Bizda etalon genom bor (pUC19 ketma-ketligi nukleotidgacha ma'lum!), shuning uchun referens bo'yicha alignment (hizalash) usulidan foydalanaman.

- Har bir tozalangan «read»ni olib, uni etalonga «shablon»dek «qo'yaman» — aligner yordamida (masalan, BWA yoki Bowtie2; ular bilan keyingi loyihalarda tanishasan).
- SAM/BAM fayllarini olaman — ular har bir «read» etalon plazmidaning aynan qayerida joylashganini va farqlar bor-yo'qligini ko'rsatadi.

Endi eng mas'uliyatli lahza — variant calling (variatsiyalarni aniqlash)! Maxsus vositalar (GATK kabi, ular bilan ham keyin tanishasan) bu fayllarni skanerlab, «read»ning har bir nukleotidini etalon bilan solishtiradi.

«Xatolar»ni qidiraman:

- SNP (snip): harfning nuqtaviy almashinuvi (A -> T)mi?
- Haqiqiy mutatsiyami yoki sekvenator xatosimi (Q-score va qoplama chuqurligiga (coverage depth) qaraymiz)?
- Yoki «xato» kattaroqmi? Unda Indel bilan ishlayapmiz: kichik insertsiya yoki deletsiya? Klonlash paytida sodir bo'lishi mumkin.

Topilgan barcha farqlar VCF hisobot fayliga yoziladi.

Variatsiyalarni filtrlayman: past sifatli hududlardagi yoki kam qo'llab-quvvatlanganlarini tashlab yuboraman. Asosiyni qidiraman: AmpR genidagi (chidamlilikni buzadimi?), ori'dagi (kopiyalar sonini pasaytiradimi?), MCS'dagi (insertni kesib olishga xalaqit beradimi?) mutatsiyalar va, albatta, MCS'dagi insertning o'zi. U kutilganga mos keladimi? Unda stop-kodonlar yoki mutatsiyalar yo'qmi?

Vaqt siqib kelmoqda…

**15:30. Haqiqat (deyarli) ochildi**

Dastlabki natijalar tayyor!

Insert MCS'da bor, uning o'lchami va ketma-ketligi asosan kutilganga mos. Insert chetida bir nechta gumonli SNP topildi — biologlar ularni qayta tekshirsin, ehtimol bu xatodir. pUC19 genlari (AmpR, lacZ fragmenti) me'yorda.

Plazmidaning GC-tarkibi ~51% — me'yor.

Qisqa hisobot tuzaman. Mukammal emas, lekin kritik muammolar yo'q! Biologlarning eksperimentini muvaffaqiyatli deb hisoblash mumkin.

**15:55. Natijalarni topshirish**

Biologlarga xulosa yuboraman: «Insert tasdiqlandi, ketma-ketlik kutilganga 98,5% mos keladi. Chetdagi ikkita nuqtaviy almashinuvni qayta tekshirish kerak. Qolgan jihatlari bo'yicha plazmida funksional. Omad!».

Pochtaga xat tushadi: «RAHMAT!!!».

Nafas rostlayman. Bugun men yana DNK-detektiv bo'lib ishladim — muddat va hamkasblar umidi bosimi ostida hayot yashirin kodini qadam-baqadam shifrladim.

Bu tergov tugadi. Oldinda — yangi ish, yangi ma'lumotlar, A, T, G va C dan iborat yangi «dalillar».

Hozircha esa — munosib tushlik.

**_Bioinformatika tibbiyot va qishloq xo'jaligini qanday o'zgartiryapti?_**

**_Genomik bioinformatikaning amaliy keyslarini Materials papkasida topasan._**

## Chapter II

Endi bioinformatika usullarini amalda qo'llash navbati senda. Sekvenirlashning real ma'lumotlarini tahlil qilish uchun olgan bilim va vositalaringdan foydalan. Ma'lumotlar sifati, adapterlar, GC-tarkibning ahamiyati va o'zingning «biologik fleshka»ngni tekshirish zarurligini unutma. Omad tilaymiz!

Ish uchun senga kerak bo'ladi:

1. O'rnatilgan Python (3.8+ versiya tavsiya etiladi).
2. O'rnatilgan [Biopython](https://biopython.org/wiki/Download) ([hujjatlar bu yerda](https://biopython.org/docs/1.75/api/index.html)).
3. [FastQC dasturi](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
4. [Trimmomatic dasturi](http://www.usadellab.org/cms/?page=trimmomatic).

### 1-topshiriq. Sekvenirlash sifatini baholash

Senga haqiqiy sekvenirlash sifatini FastQC dasturi yordamida baholash kerak bo'ladi.

1. Papkalar strukturasini yarat:

- project1
  - results
  - data
  - script

2. [NCBI](https://www.ncbi.nlm.nih.gov/sra/SRX1684458%5Baccn%5D) sahifasini och. Bu yerda sekvenirlash haqidagi ma'lumotni ko'rasan. Ma'lumotni o'rgan. Keyin sahifadan Run identifikatorini top. Unga bosib, assembly sahifasiga o'tasan. Senga `.fastq` kengaytmali faylni yuklab olish kerak. Uni data papkasiga joylashtir.
3. FastQC yordamida sekvenirlash ma'lumotlarini tahlil qil. Barcha olingan natija fayllarini results papkasiga joylashtir, fayl nomini topshiriq band va kichik bandiga muvofiq yoz (masalan, `1.1.txt`).
4. FastQC tahlilida qanday kritik muammolar (qizil status) aniqlay olganingni belgila:

   - Basic Statistics.
   - Per base sequence quality.
   - Per tile sequence quality.
   - Per sequence quality scores.
   - Per base sequence content.
   - Per sequence GC content.
   - Per base N content.
   - Sequence length distribution.
   - Sequence duplication levels.
   - Overrepresented sequence.
   - Adapter content.

   Qaysi muammolar e'tibor talab qilishini (to'q sariq status) belgila:

    - Basic Statistics.
    - Per base sequence quality.
    - Per tile sequence quality.
    - Per sequence quality scores.
    - Per base sequence content.
    - Per sequence GC content.
    - Per base N content.
    - Sequence length distribution.
    - Sequence duplication levels.
    - Overrepresented sequence.
    - Adapter content.

5. Adapterlar bo'yicha tozalikka (tegishli bo'limdagi grafik) va Overrepresented sequences hisobotiga e'tibor ber. Hisobot va grafikni tahlil qilishda topilgan adapterlar nomlarini yoz.

### 2-topshiriq. Sekvenirlash sifatini yaxshilash (reads'larni tozalash)

Endi sening vazifang — read'larni adapterlardan tozalash. Buning uchun Trimmomatic dasturidan foydalan. Barcha olingan natija fayllarini results papkasiga joylashtir, fayl nomini topshiriq raqami va kichik raqamiga muvofiq yoz.

1. Xom ma'lumotlarni adapterlardan tozalash uchun avval ularni aniqlash kerak.

   Hisobotga diqqat bilan qara. Trimmomatic dasturida tozalash uchun standart adapter fayllari mavjud (Trimmomatic o'rnatilgandan so'ng adapters papkasidan mosini top va data papkasiga joylashtir). Bu bosqichda senga ILLUMINACLIP funksiyasi kerak bo'ladi, parametrlarni mustaqil tanlash kerak — asl fayl, sekvenirlash tavsifi va [qo'llanma (manual)]ga (http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) tayanib.

2. Trimmomatic dasturi yordamida o'zing tanlagan adapterlar faylida read'lar sifatini yaxshilashga harakat qil, hosil bo'lgan faylni results papkasiga saqla.

3. FastQC'ni ishga tushir va adapterlardan tozalangan to'plamni tahlil qil, hosil bo'lgan faylni results papkasiga saqla.

4. Adapterlar tozalandimi? Agar adapterlar tozalanmagan bo'lsa, boshqa fayldan foydalanib ko'r. Eng mos adapterlar faylini data papkasiga nusxala.

5. FastQC'ni ishga tushir va adapterlardan tozalangan to'plamni tahlil qil. Xom ma'lumotlarni yakuniy tayyorlash uchun qanday kritik muammolar qolganini ko'r. Ularni belgila.
    - Basic Statistics.
    - Per base sequence quality.
    - Per tile sequence quality.
    - Per sequence quality scores.
    - Per base sequence content.
    - Per sequence GC content.
    - Per base N content.
    - Sequence length distribution.
    - Sequence duplication levels.
    - Overrepresented sequence.
    - Adapter content.

6. Tim-lidingni so'raydigan savol:  
  _«Nima uchun read'lar uzunligi taqsimoti o'zgardi?»_  
  P2P tekshiruv uchun batafsil javob tayyorla.

7. Yana asl fayl bilan Trimmomatic dasturida ishla (bu safar Trimmomatic'ning barcha mumkin bo'lgan funksiyalari bilan tanish, standart parametrlarni qo'lla, natija fayllarini (Trimmomatic, FastQC) results papkasiga saqla).

8. Tim-lidingni so'raydigan savol:  
  _«O'zing tanlagan parametrlarni tasvirla. Nima uchun aynan shunday tanlov qilindi? Hisobotni yaxshilash muvaffaq bo'ldimi (ha/yo'q) va nima uchun?»_  
  P2P tekshiruv uchun batafsil javob tayyorla.

### 3-topshiriq. GC-tarkibni baholash

Endi pUC19'ning ikki xil vektorida GC-tarkibni hisoblab ko'raylik.

1. fasta ketma-ketliklarni yuklab ol va ularni [fasta1](https://www.ncbi.nlm.nih.gov/nuccore/6691170>) va [fasta2](https://www.ncbi.nlm.nih.gov/nuccore/MT856194.1>) deb nomla.

2. Ularni data papkasiga joylashtir.

3. GC-tarkibni % da aniqlaydigan skript yoz. Tayyor skriptni script papkasiga qo'y.

4. Tim-lidingni so'raydigan savol:  
  _«pUC19 vektorining etalon GC-tarkibi qancha? Olingan qiymatlar bilan solishtir. Sening fikringcha, farq bormi va nima uchun?»_  
  P2P tekshiruv uchun batafsil javob tayyorla.

💡 [Bu yerni bos](https://new.oprosso.net/p/4cb31ec3f47a4596bc758ea1861fb624), **shu loyiha bo'yicha bizga fikr-mulohaza bildirish uchun**. So'rovnoma anonim va jamoamizga o'qitishni yaxshilashda yordam beradi. Loyihani tugatgan zahoti to'ldirishni tavsiya qilamiz.
