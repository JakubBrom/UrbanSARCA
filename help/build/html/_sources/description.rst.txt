Description of Urban Green SARCA (In Czech)
============================================

Modul Urban Green SARCA je softwarový GIS nástroj vytvořený pro účely odhadu
depozice radionuklidu na povrchu vegetace a půdy, respektive ploch bez
vegetace v časné fázi radiační havárie. Modul umožňuje odhadnout kontaminaci
vegetačního krytu prostřednictvím satelitních snímků, podkladů o celkové
depozici radionuklidu a na základě informace o úhrnu srážek v průběhu
depozice radionuklidu v zájmovém území. Výpočet je tedy prováděn pro podmínky
suché i mokré depozice radionuklidu. Vzhledem k tomu, že modul umožňuje
výpočet jednotlivých proměnných s využitím podrobných satelitních dat, je
využití modulu vhodné zejména pro urbánní území a pro území, kde nejsou k
dispozici informace o jednotlivých typech zeleně (vegetace) nebo je
identifikace zeleně značně komplikovaná.

Modul je koncipován tak, aby minimalizoval množství vstupů a zároveň
poskytoval dostatečné množství výstupů důležitých pro následné rozhodování v
oblasti radiační ochrany urbánních území.


Funkcionalita modulu Urban Green SARCA
---------------------------------------

Funkcionalita modulu vychází z programu SARCA (Brom et al. 2015), který je
nicméně určen pro analýzu kontaminace polních plodin a je založen na
modelování časových změn produkčních charakteristik plodin. Oproti tomu,
Urban Green SARCA využívá pro odhad produkčních charakteristik vegetace
dostupná satelitní data.

Výpočet jednotlivých výstupů lze rozdělit do dvou částí. V první části jsou
ze satelitních dat vypočteny produkční charakteristiky, tedy množství biomasy
a index listové plochy vegetace. Ve druhé části je vypočtena depozice
radioaktivního kontaminantu na povrchu vegetace a na ostatních površích. Na
základě vypočtených dat jsou dále hodnoceny kategorie radioaktivní
kontaminace vegetace (zeleně) podle stanovených referenčních úrovní
kontaminace (viz dále). Doplňkově je vypočtena vrstva intercepčního faktoru,
tedy relativní distribuce radionuklidu mezi zeleň a ostatní povrchy a
vypočtena je též hmotnostní kontaminace, včetně vyznačení nadlimitní
hmotnostní kontaminace zeleně. Způsoby výpočtu jednotlivých ukazatelů jsou
uvedeny dále.


Výpočet výstupů modulu Urban Green SARCA
-----------------------------------------

Množství biomasy
.................

Při výpočtu množství biomasy (B; :math:`t.ha^{-1}`) na dané ploše vycházíme ze
vztahu mezi množstvím zelené biomasy a jejím spektrálním projevem. Pro
vyjádření byl použit následující zjednodušený vztah pro odhad množství biomasy:

.. math::

    B=50\cdot NDVI^{2.5},


kde NDVI je normalizovaný rozdílový vegetační index (Normalized Difference
Vegetation Index; Rouse Jr et al., (1973)):

.. math::

    NDVI=\frac{R_{NIR}-R_{RED}}{R_{NIR}+R_{RED}},

kde :math:`R_{NIR}` a :math:`R_{RED}` jsou spektrální reflektance v blízké
inrfačervené (NIR; přibližně 800 nm) a v červené oblasti (RED; přibližně 670
nm) (rel.). Použitý model představuje hrubý odhad živé biomasy zeleně. Do
budoucna předpokládáme nahrazení uvedeného vztahu komplexnějším modelem.


Index listové plochy
.....................

Index listové plochy (bezrozm.) je počítán pomocí jednoduchého lineárního
vztahu mezi listovou plochou a spektrálním indexem NDVI:

.. math::

    LAI=4.9\cdot NDVI-0.46

Kontaminace zeleně a půdy, intercepční faktor
..............................................

Pro rozhodování o množství depozice radioaktivního materiálu na povrchu
porostu a povrchu půdy je vypočten intercepční faktor (rel.), který je
ukazatelem, jak velká frakce depozice zůstává na povrchu porostu. Hodnota
závisí na indexu listové plochy porostu a úhrnu srážek v průběhu depozice.
Podle Müllera a Pröhla (1993) lze intercepční frakci (faktor) depozice
radioizotopu fw v časné fázi radiační havárie vypočítat podle vzorce:

.. math::

    f_{w}=\min\left[1;\frac{LAI\cdot k\cdot S\left(1-\mathrm{e^{-\frac{\ln2}{3S}R}}\right)}{R}\right]

kde k je specifický faktor pro daný kontaminant (I: k = 0.5; Sr, Ba: k = 2;
Cs a ostatní radionuklidy: k = 1), S je tloušťka vodního filmu na rostlinách
(mm) a R je úhrn srážek (mm). Hodnota S je zpravidla 0,15 – 0,3 mm se střední
hodnotou 0,2 mm (Pröhl, 2003). Výpočet depozice na povrchu rostlin vychází z
předpokladu, že depozice na povrchu rostlin je poměrnou částí celkové
depozice danou intercepčním faktorem:

.. math::

    D_{biomasa}=D_{celk}\cdot f_{w}

kde :math:`D_{biomasa}` je měrná depozice radioizotopu na povrchu rostlin
:math:`(Bq.m^{-2})` a Dcelk je celková měrná radioaktivní depozice :math:`(Bq
.m^{-2})` zadávaná jako vstup do modelu. Měrná depozice radioizotopu na
povrchu půdy (Dpuda ; :math:`Bq.m^{-2}`) je pak rozdílem mezi celkovou měrnou
depozicí a měrnou depozicí na povrchu porostu:

.. math::

    D_{puda}=D_{celk}-D_{biomasa}

Pokud jsou hodnoty vypočteného množství biomasy menší než 0,5 :math:`t
.ha^{-1}`, je vypočtena pouze měrná depozice radioaktivního materiálu na
povrchu půdy. Důvodem je minimální předpoklad možnosti odstranění biomasy.
Doplňkovou charakteristikou je hmotnostní kontaminace biomasy zeleně (Dhmot;
:math:`Bq.m^{-2}`), která je vypočtena podle vztahu:

.. math::

    D_{hmot}=\frac{D_{biomasa}}{B \cdot 0.1}


Referenční úrovně a  překročení hygienického limitu kontaminace biomasy
.........................................................................

Území kontaminované radioaktivní depozicí je pro praktické účely rozděleno na
tři oblasti, v závislosti na stanovených referenčních úrovních. Rozdělení
sledovaného území do oblastí podle referenčních úrovní vychází z předpokladu,
že lze vymezit území, ve kterých kontaminace nepřekračuje stanovenou úroveň
dávkového příkonu nebezpečného pro obyvatelstvo a zvířata (hodnota 0), dále
území ve kterých lze provádět opatření za účelem radiační ochrany (hodnota 1)
a území, kde úroveň radioaktivní kontaminace, respektive dávkového příkonu
překračuje bezpečnou hranici pro další management (hodnota 2).
Pro referenční úrovně RU 0 a RU 2 není doporučeno odstranění biomasy za
účelem ochrany půdy. V prvním případě (RU 0) nepřesahuje kontaminace
stanovenou mez a nejsou ze předpokládána další rizika, zeleň a produkci
rostlinné biomasy je možné využít běžným způsobem, případně v omezené míře na
základě dalších postupů. Naopak v případě ploch zařazených do referenční
úrovně RU 2 existuje předpoklad nadlimitní radioaktivní kontaminace ploch a
možnost ohrožení zdraví pracovníků pověřených manipulací s nadzemní biomasou
rostlin. V rámci ploch zařazených do RU 1 lze předpokládat půdoochranný
význam vegetačního krytu, který lze za daných podmínek odstranit z půdního
povrchu. Limitem je zde množství živé nadzemní biomasy 0,5 :math:`t.ha^{-1}`,
kdy předpokládáme, že sklizeň menšího množství biomasy na danou plochu je již
neefektivní, případně technicky nemožná. Hranice referenčních úrovní lze
nastavit přímo v uživatelském rozhraní Urban Green SARCA. Přednastaveny jsou
hodnoty 5000 :math:`Bq.m^{-2}` a 3 :math:`MBq.m^{-2}`.
Vedle vrstvy referenčních úrovní je výstupem modelu vrstva překročení
hygienického limitu kontaminace biomasy. Rastrová vrstva nese hodnoty 0 pro
pixely, ve kterých je zjištěna hodnota úrovně hmotnostní kontaminace biomasy
menší než stanovená hodnota v uživatelském rozhraní Urban Green SARCA.
Hodnoty přesahující stanovený hygienický limit jsou zařazeny do kategorie 1.
V modulu Urban Green SARCA je přednastavena hodnota 1000 :math:`(Bq.kg^{-1})`,
která odpovídá nejvyšší přípustné úrovni radioaktivní kontaminace potravin pro
radiačně mimořádné situace podle vyhlášky 389/2010 Sb. o radiační ochraně
(Vyhláška 389/2012 Sb. o radiační ochraně, 2012).


Přehled použité literatury
--------------------------

*Muller, H., Prohl, G., 1993. Ecosys-87: A dynamic model for assessing
radiological consequences of nuclear accidents. Health Phys. 64, 232–252.*

*Pröhl, G., 2003. Radioactivity in the terestrial environment, in: Scott, E.M.
(Ed.), Modelling Radioactivity in the Environment. Elsevier, Amsterdam;
Boston, pp. 87–108.*

*Rouse Jr, J., Haas, R., Schell, J., Deering, D., 1973. Monitoring vegetation
systems in the Great Plains with ERTS In Third Earth Resources Technology
Satellite-1, in: Third Earth Resources Technology Satellite-1 Symposium: The
Proceedings of a Symposium Held by Goddard Space Flight Center at Washington,
DC on December 10-14, 1973: Prepared at Goddard Space Flight Center.
Scientific and Technical Information Office, NASA, pp. 309–317.*

*Vyhláška 389/2012 Sb. o radiační ochraně, 2012.*