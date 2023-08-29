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

Index listové plochy (bezrozm.) je počítán na základě
spektrálních dat v červené a blízké infračervené oblasti. Jedná se o
jednoduchý přístup, který nicméně dobře koresponduje především
indexem listové plochy bylinné vegetace, zejm. polních plodin.
Do modelu byla přidána možnost volby metody výpočtu indexu
listové plochy. K dispozici jsou následující metody:

**1. Jednoduchá metoda (Simple method)**

Metoda využívá lineární vztah mezi listovou plochou a spektrálním
indexem NDVI:

.. math::

    LAI=4.9 \cdot NDVI-0.46

Metoda poskytuje vcelku univerzální výsledky pro nízkou i vysokou
lesní vegetaci. Ostatní jmenované metody jsou vhodné spíše pro polní
kultury

**2. Pôças**

Podle Pôças et al. (2014).

.. math::

    LAI =
      \begin{cases}
        11 \cdot SAVI^3         & SAVI > 0;\  SAVI \leq 0.817\\
        6                       & SAVI > 0.817
      \end{cases}

kde :math:`SAVI` je spektrální index (Soil Adjusted Vegetation Index;
Huete, 1988):

.. math::
    SAVI = \frac{(1 + L) \cdot (R_{NIR} - R_{RED})}{L + R_{NIR}+R_{RED}}

kde :math:`L` je konstanta (:math:`L=0.5`).


**3. Bastiaanssen**

Podle Bastiaanssena et al. (1998)

.. math::

    LAI =
      \begin{cases}
        -\frac{\ln \frac{0.61-SAVI}{0.51}}{0.91}
            & SAVI >0;\  SAVI\leq 0.61\\
        6                       & SAVI > 0.61
      \end{cases}

**4. Jafaar**

Podle Jafaar et al. (2019).

.. math::
    LAI_1=
        \begin{cases}
            11 \cdot SAVI^3         & SAVI > 0;\  SAVI \leq 0.817\\
            6                       & SAVI > 0.817
        \end{cases}

.. math::
    LAI_2 =
      \begin{cases}
        -\frac{\ln \frac{0.61-SAVI}{0.51}}{0.91}
            & SAVI >0;\  SAVI\leq 0.61\\
        6                       & SAVI > 0.61
      \end{cases}

.. math::
    LAI = \frac{LAI_1 + LAI_2}{2}

**5. Anderson**

Podle Andersson et al. (2004).

.. math::
    LAI = (4\cdot OSAVI -0.8) \cdot (1 + 4.73 \cdot 10^{-6}\mathrm{e}^{15.64 \cdot OSAVI})

kde :math:`OSAVI` je spektrální index (Optimized Soil Adjusted
Vegetation Index; Rondeaux et al., 1996):

.. math::
    OSAVI = \frac{R_{NIR} - R_{RED}}{0.16 + R_{NIR}+R_{RED}}

**6. Carrasco**

Podle Carrasco-Benavides et al. (2014).

.. math::
    LAI = 1.2 - 3.08\mathrm{e}^{-2013.35 \cdot NDVI^{6.41}}

**7. Turner**

Podle Turner et al. (1999).

.. math::
    LAI=0.5724 + 0.0989 \cdot NDVI - 0.0114 \cdot NDVI^2 + 0.0004\cdot NDVI^3


**8. Haboudane**

.. math::
    LAI = 0.0918^{6.0002 \cdot RDVI}

kde :math:`RDVI` je spektrální index (Renormalized Difference
Vegetation Index; Roujean and Breon, 1995):

.. math::
    RDVI = \frac{R_{NIR} - R_{RED}}{\sqrt{R_{NIR} + R_{RED}}}

**9. Brom**

.. math::
    LAI = \frac{6.0}{1 + \mathrm{e}^{-(8 \cdot SAVI - 5)}}

Uvedený přístup je experimentální. Metoda velmi dobře koresponduje s
přístupy podle Bastiaanssena, Jafaara a Pôçase, nicméně je velmi
senzitivní na hodnoty indexů v exponentu funkce.


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
nastavit přímo v uživatelském rozhraní Urban Green SARCA.
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

*Anderson, M., Neale, C., Li, F., Norman, J., Kustas, W., Jayanthi,
H., Chavez, J., 2004. Upscaling ground observations of vegetation
water content, canopy height, and leaf area index during SMEX02 using
aircraft and Landsat imagery. Remote Sensing of Environment 92,
447–464. https://doi.org/10.1016/j.rse.2004.03.019*

*Bastiaanssen, W.G.M., Menenti, M., Feddes, R.A., Holtslag, A.A.M.,
1998. A remote sensing surface energy balance algorithm for land (
SEBAL). 1. Formulation. Journal of Hydrology 212–213, 198–212.
https://doi.org/10.1016/S0022-1694(98)00253-4*

*Carrasco-Benavides, M., Ortega-Farías, S., Lagos, L., Kleissl, J.,
Morales-Salinas, L., Kilic, A., 2014. Parameterization of the
Satellite-Based Model (METRIC) for the Estimation of Instantaneous
Surface Energy Balance Components over a Drip-Irrigated Vineyard.
Remote Sensing 6, 11342–11371. https://doi.org/10.3390/rs61111342*

*Huete A.R. (1988): A soil-adjusted vegetation index (SAVI) Remote
Sensing of Environment 27, 47-57.*

*Jaafar, H.H., Ahmad, F.A., 2019. Time series trends of Landsat-based
ET using automated calibration in METRIC and SEBAL: The Bekaa
Valley, Lebanon. Remote Sensing of Environment S0034425718305947.
https://doi.org/10.1016/j.rse.2018.12.033*

*Muller, H., Prohl, G., 1993. Ecosys-87: A dynamic model for assessing
radiological consequences of nuclear accidents. Health Phys. 64, 232–252.*

*Pôças, I., Paço, T.A., Cunha, M., Andrade, J.A., Silvestre, J.,
Sousa, A., Santos, F.L., Pereira, L.S., Allen, R.G., 2014.
Satellite-based evapotranspiration of a super-intensive olive
orchard:  Application of METRIC algorithms. Biosystems Engineering
128, 69–81. https://doi.org/10.1016/j.biosystemseng.2014.06.019*

*Pröhl, G., 2003. Radioactivity in the terestrial environment, in: Scott, E.M.
(Ed.), Modelling Radioactivity in the Environment. Elsevier, Amsterdam;
Boston, pp. 87–108.*

*Rondeaux G., Steven M., Baret F. (1996): Optimisation of
soil-adjusted vegetation indices Remote Sensing of Environment,
55 (1996), pp. 95-107*

*Roujean, J.-L., Breon, F.-M., 1995. Estimating PAR absorbed
by vegetation from bidirectional reflectance measurements.
Remote Sensing of Environment 51, 375–384.
https://doi.org/10.1016/0034-4257(94)00114-3*

*Rouse Jr, J., Haas, R., Schell, J., Deering, D., 1973. Monitoring vegetation
systems in the Great Plains with ERTS In Third Earth Resources Technology
Satellite-1, in: Third Earth Resources Technology Satellite-1 Symposium: The
Proceedings of a Symposium Held by Goddard Space Flight Center at Washington,
DC on December 10-14, 1973: Prepared at Goddard Space Flight Center.
Scientific and Technical Information Office, NASA, pp. 309–317.*

*Turner, D.P., Cohen, W.B., Kennedy, R.E., Fassnacht, K.S., Briggs,
J.M., 1999. Relationships between Leaf Area Index and Landsat TM
Spectral Vegetation Indices across Three Temperate Zone Sites.
Remote Sensing of Environment 70, 52–68.
https://doi.org/10.1016/S0034-4257(99)00057-7*

*Vyhláška 389/2012 Sb. o radiační ochraně, 2012.*