#import "@preview/charged-ieee:0.1.3": ieee

#show: ieee.with(
  title: [Proteomics],
  abstract: [
    In dit project gaan we RAW-data van een reeks massaspectrometrische metingen verwerken. Hiervoor wordt MSfragger gebruikt om een omzetting naar pepXML uit te voeren. Om peptide spectrum matches (PSM’s) binnen deze pepXML bestanden te analyseren, zijn er verschillende methoden ontwikkeld in python en R. 
  ],
  authors: (
(

      name: "Anass Hamzaoui",
      department: [Faculty of Science],
      organization: [University of Antwerp],
      location: [Antwerp, Belgium],
      email: "Anass.hamzaoui@student.uantwerpen.be",
    ),
    (
      name: "Chan Min Jan",
      department: [Faculty of Science],
      organization: [University of Antwerp],
      location: [Antwerp, Belgium],
      email: "Chan.min.jan@student.uantwerpen.be",
    ),
(
      name: "Axel De Leeuw",
      department: [Faculty of Science],
      organization: [University of Antwerp],
      location: [Antwerp, Belgium],
      email: "Axel.de.leeuw@student.uantwerpen.be",
    ),
    
    
),
  index-terms: ("Scientific writing", "Typesetting", "Document creation", "Syntax"),
  bibliography: bibliography("refs.bib"),
  figure-supplement: [Fig.],
)

= Introduction

In ons onderzoeks project hebben wij ons gericht op het veld van proteomics/ spectrometrie, waaronder post-translationele modificaties (PTMs) en aminozuursubstituties. Deze variaties beïnvloeden eiwitfuncties en kunnen verband houden met gezondheid en ziekte. Dit onderzoek richt zich op twee hoofdvragen: (1) Of een Open Search-benadering meer spectra identificeert dan een traditionele Closed Search, en (2) of PTM-patronen en aminozuursubstituties verschillen tussen etnische groepen, zoals White versus Black/African American populaties.

Met FragPipe (inclusief MSFragger) verwerkten we MS-data tot geannoteerde eiwitprofielen, waarna Python en R werden gebruikt voor statistische en bio-informatica-analyses. Door spectra-annotatiemethoden te vergelijken, willen we inzicht krijgen in de voordelen van Open Search. Daarnaast onderzoeken we of genetische of omgevingsfactoren leiden tot verschillen in PTMs tussen populaties.

Alle gebruikte code en analysemethoden zijn openbaar beschikbaar in onze [#link("https://github.com/Axelmans/Proteomics")[Github Repository]], zodat onze bevindingen reproduceerbaar en transparant zijn. Dit onderzoek draagt bij aan een beter begrip van eiwitvariatie en kan toekomstige studies naar genetische en moleculaire mechanismen ondersteunen.

= Dataset
Voor dit onderzoek werd een selectie gemaakt van tumorstalen afkomstig uit verschillende etnische groepen: African, Asian, Hispanic, Native, Other, en White. Deze dataset biedt een unieke mogelijkheid om proteomische variaties te bestuderen in relatie tot genetische achtergrond, wat relevant is voor onderzoek naar gepersonaliseerde geneeskunde en tumorbiologie.

De data bestaat uit .mzML-bestanden (geconverteerde massaspectrometrie-rawfiles) per etnische groep, samen met een decoy-database in de vorm van een fasta.fas-bestand. Deze bestanden vormen de basis voor eiwitidentificatie en -kwantificatie met behulp van FragPipe (MSFragger, Philosopher, en andere tools). Het gebruik van een decoy-database helpt bij het minimaliseren van false discovery rates (FDR), wat essentieel is voor betrouwbare resultaten in grootschalige proteomics-analyses.

Een belangrijk aspect van deze dataset is de etnische diversiteit, waardoor we verschillen in post-translationele modificaties (PTMs) en aminozuursubstituties tussen populaties kunnen onderzoeken. Kleine genetische variaties tussen etnische groepen kunnen leiden tot verschillen in eiwitexpressie of -structuur, wat mogelijk invloed heeft op ziekteprogressie of therapierespons. Door deze dataset te heranalyseren, streven we ernaar om nieuwe inzichten te genereren die kunnen bijdragen aan precisiegeneeskunde, waarbij behandelingen beter afgestemd kunnen worden op individuele (en populatie-specifieke) kenmerken.

De combinatie van massaspectrometrie-gegevens en bio-informatica-analyses (Python/R) stelt ons in staat om zowel globale trends als subtiele, maar klinisch relevante, verschillen tussen populaties te detecteren. Dit maakt de dataset niet alleen waardevol voor fundamenteel onderzoek, maar ook voor toekomstige translationele toepassingen.


#v(10pt)
= Closed vs. Open search
#v(10pt)
Closed Search en Open Search zijn beide technieken om spectra te matchen aan peptiden. Het voornaamste verschil is dat bij Closed Search er enkel gematched wordt indien het massaverschil tussen het spectrum en peptide zeer klein is (b.v. ±20 ppm).
Dit heeft als gevolg dat enkel spectra met weinig tot geen modificaties gematched worden.
Bij Open Search mag dit verschil veel groter zijn, wat veel meer (combinaties van) modificaties toestaat op spectra terwijl deze nog steeds gematched kunnen worden.

#figure(
  image("ReportData/Picture26.png")
)


Het is echter niet gegarandeerd dat Open Search hierdoor meer spectra annoteert.                         In het geval dat alle spectra een klein massaverschil vertonen met een bepaald peptide draagt Open Search niets bij, terwijl dat voor dit onderzoek zeker wel gewenst is.                    Daarom werd vooraf een vergelijkend onderzoek gedaan tussen de 2 op onze dataset.

De code gaat elk spectrum af, en houdt enkel rekening met betrouwbare matches.              Hiervoor werd gekeken naar de zogenaamde ‘expect’ score, die de kans geeft dat een match louter toeval was, bij voorkeur is deze kans laag (<= 5%). Tussen alle betrouwbare matches werd onderscheid gemaakt tussen target matches en decoy matches om de False Discovery Rate (FDR) te kunnen bepalen, ook deze waarde is bij voorkeur klein (<= 5%).

De vergelijking tussen de 2 technieken gaf de volgende conclusies:

1.	Bij Open Search neemt het aantal betrouwbare annotaties toe met 50% ten opzichte van het aantal annotaties bij Closed Search (1838696 vs. 1253250).
2.	Het aantal decoy matches is iets groter bij Closed Search (10761 vs. 9820).
3.	Closed Search heeft een grotere FDR, maar blijft onder 5% (0,009 vs. 0,005).

Een Open Search uitvoeren is in dit geval dus zinvol, zoals gewenst.

#v(10pt)
== Mogelijke uitbreiding en verbetering
#v(10pt)

Dit onderzoek kan uiteraard nog verbeterd en/of uitgebreid worden:

1.	De vergelijking tussen Closed en Open Search zou uitgebreid kunnen worden:                 het zou bijvoorbeeld interessant kunnen zijn om de overlap in peptiden die gematched werden te vergelijken tussen de 2 technieken.
2.	Naast PTM’s te kunnen matchen zou het ook interessant zijn om te kunnen bepalen waar deze PTM’s zich exact op het peptide bevinden, dit viel echter niet te achterhalen uit de data in de .pepXML bestanden.
3.	De opschaling in Python kon eventueel afgerond worden met statistische significantietesten die de waarnemingen bevestigen of ontkrachten.
4.	Om de etnische groepen verder te vergelijken, zou het ook interessant kunnen zijn om aminozuursubstituties te achterhalen en hoe dit mogelijks verschilt per groep.



= PTM analyze
== R implementatie
Het doel van dit R script is een om een diagnostische tool te ontwikkelen om specifieke regulatorische pathways te ontdekken op basis van post-translationele modificaties (PTM’s). Dit script bevat een uitgebreide lijst aan gekende PTM’s die in diverse fysiologische condities relevant kunnen zijn.
Link naar het script: #link("https://github.com/Axelmans/Proteomics/tree/main/rcode")[R script]
#v(10pt)
=== Core Logic
#v(10pt)
Eerst wordt er, via een for-loop, uit de bemonsterde PSM’s van de pepXML-file een verschil berekend tussen de precursor molecule en het herkende peptide fragment, dit verschil wordt vervolgens opgeslagen in een vector: 

#figure(
  image("ReportData/CoreLogic1.png")
)
#figure(
  image("ReportData/CoreLogic2.png")
)

#figure(
  image("ReportData/CoreLogic3.png")
)

Vervolgens wordt de PTM_matcher functie opgeroepen om met een zo goed mogelijke PTM combinatie, Δm te benaderen. Er wordt dus eerst een reeks combinaties voorgesteld die dan beoordeeld worden of ze al dan niet een goede fit zijn voor Δm. Het genereren van PTM combinaties wordt gedaan vanuit de PTM lijst impliciet in de matcher functie, het gaat als volgt:

#figure(
  image("ReportData/Picture4.png")
)

Dus als er bijvoorbeeld in de PTM lijst slechts drie PTM’s zijn: A, B, en C 
dan worden er alle mogelijke combinaties gegenereerd met combn als volgt: 

#figure(
  image("ReportData/Picture5.png")
)

Hun overeenkomstige massas worden dan ook berekend, deze zijn voorbepaald in de lijst met PTM’s. Na berekening worden ze opgeslagen in de variabele combo_masses en worden ze benoemd met names (dit is handig voor de output). Vervolgens gaan we de beste fit/match bepalen door te kijken welke combinatie het kleinste overschot maakt. Dit wordt nagegaan door een verschil te maken tussen de gekozen PTM-combinatie en de Δm, het getal wordt dan beoordeeld ten opzichte van een tolerantie (έ): 

#figure(
  image("ReportData/Picture6.png", width: 60%)
)

#figure(
  image("ReportData/Picture7.png", width: 80%)
)

Wanneer echter de gefitte som de tolerantie overschrijdt, word er een nieuwe poging tot geschikte combinatie gemaakt (trial) met een verdubbelde tolerantie. De matches kunnen dan gerankt worden naargelang hun aantal trials, hoe meer trials er werden uitgevoerd hoe minder betrouwbaar de match. 

Na 5 pogingen de stopt de code want de PTM combinaties die na 5 maal de tolerantie gefit worden zullen niet betrouwbaar zijn, het verschil tussen de potentiële match en Δm is te groot om betekenisvol te zijn. diffs is het verschil tussen de fit en Δm, deze zal in de output delta_fit genoemd worden. length(best_idx) > 0 wilt zeggen dat er een match is gevonden, als er geen match is gevonden is length(best_idx) == 0 en dan krijgen we een NA als output. Met best < - best_idx[which.min(diffs[best_idx])] wordt degene met de kleinste diffs gekozen als match en word er een dataframe gegeven als output.  
#v(10pt)
=== pepXML subsampling strategy
#v(10pt)
Alle PSM’s worden opgezocht en de totale hoeveelheid daarvan wordt geteld. Gezien pepXML files enorm veel PSM’s bevatten werd er gekozen om een bemonstering te doen. Op 5 verschillende plaatsen in de totale hoeveelheid aan PSM’s (chunks) worden er 150 stuks genomen. Deze 5 plaatsen zijn gelijkmatig verdeeld door simpelweg de totale hoeveelheid PSM’s te delen door 5.

De kans dat bijvoorbeeld carbamylatie en methylatie onafhankelijk voorkomen is bijvoorbeeld extreem klein. De kleinste p-waarden heb ik van boven laten verschijnen met order. 
#figure(
  image("ReportData/Picture8.png")
)

De redenering is dat er dan een zo representatief mogelijke bemonstering wordt gedaan dankzij de gelijkmatige verdeling. Dit werd voornamelijk gedaan uit gebrek aan rekenkracht. Hoeveel chunks je wil nemen en hoeveel PSM’s je daaruit wilt halen kan makkelijk worden veranderd. In de code: sample(start_idx:end_idx, sample_size) wordt er gezorgd dat er binnen de chunk random wordt gekozen. sample_indices kan je zo nodig uitprinten om te zien welke PSM’s je uit de pepXML hebt gehaald, deze variabele bevat alle indices van de PSM’s die gesampled worden binnen de chunks.
Onderaan is bijvoorbeeld een print van de sample_indices bij een subsampling van 5 chunks waaruit 25 random PSM’s worden geselecteerd. De pepXML file bevat 32700 PSM’s. 

#figure(
  image("ReportData/Picture9.png")
)

#figure(
  image("ReportData/Picture10.png", height: 40%)
)
#v(10pt)
=== Output example
#v(10pt)
overlopen van de kolommen : 
-	peptide : fragment-ionen die gereconstrueerd zijn tot een peptide-sequentie 
-	protein : het eiwit van waar de sequentie afkomstig is 
-	precursor_mass : het moleculair gewicht van de precursor-ion vooraleer het doorheen de collisiekamer ging 
-	pepfrag_mass : massa van de herkende peptide-sequentie
-	delta_mass : het verschil in moleculair gewicht tussen de herkende peptide-sequentie uit de fragmentionen en de oorspronkelijke moedermolecule (Δm). 
-	PTM_combo_mass : het moleculair gewicht van de gefitte som 
-	delta_fit : het verschil tussen de gefitte som en pepfrag_mass (diffs in de code)
-	match : de keuze aan PTM’s die gemaakt werden om in delta_mass te passen 
-	used_tolerance  : de gebruikte tolerantie, het maximaal toegelaten verschil (έ) tussen PTM_combo_mass en pepfrag_mass 
-	trials : hoeveel keer is de fit van de PTM_combo gefaald of hoe vaak is de tolerantie vergroot moeten worden om een gepaste combinatie te zoeken. Hoe meer trials er zijn gedaan hoe minder betrouwbaar de fit
#v(10pt)
=== Diagnostic plots
#v(10pt)
We kunnen het voorkomen van bepaalde PTM’s binnen onze bemonstering bekijken. Hieruit kunnen eventueel initiële aanwijzingen voor fysiologische condities gevonden worden (bv. Immuunrespons of signaaltransducties voor differentiële genexpressies).
Onderstaande plot weergeeft de frequentie van afzonderlijke PTM’s per spectrum querie.

#figure(
  image("ReportData/Picture11.jpg")
)

Anderzijds kunnen we ook het voorkomen van de massas van de PTM combinaties bekijken voor gelijkaardige aanwijzingen. 

#figure(
  image("ReportData/Picture12.jpg")
)
Om een beeld te hebben van hoe adequaat de PTM’s zijn gefit zijn ten opzichte van Δm  zijn de frequenties van gebruikte toleranties  uitgeplot, zie onderaan. 

#figure(
  image("ReportData/Picture13.jpg")
)
#v(10pt)
=== Co-occurence analysis with multidimensional scaling
#v(10pt)
Hier gaan we kijken welke PTM’s vaak samen voorkomen. Door te kijken welke PTM’s clusters vormen kunnen we aannemen dat ze deel uitmaken van een signaaltransductie, hieruit kunnen dan beweringen worden gemaakt met betrekking tot de fysiologische staat van het oorspronkelijk biotisch staal. 

De match strings van ons resultaat final_result gaan we opdelen in individuele PTM’s en we halen hieruit alle PTM’s zonder duplicaten met unique. PSM’s zonder match gaan we negeren met na omit. 

#figure(
  image("ReportData/Picture14.png")
)

Nu hebben we allemaal afzonderlijke PTM’s waarvan we een binaire matrix van kunnen opstellen, dit is nodig om een Jaccard afstand te bepalen 

#figure(
  image("ReportData/Picture15.png")
)

DIt creërt een matrix zoals : 

#figure(
  image("ReportData/table1.png")
)

In R is 0 en 1 TRUE of FALSE. Deze waarden kunnen worden gebruikt om de Jaccard afstand te berekenen. Jaccard similariteit is de verhouding tussen de doorsnee en unie van 2 verzamelingen. De Jaccard afstand binnen onze context is uitgedrukt als volgt: 

#figure(
  image("ReportData/Picture16.png")
)

Dit wordt pas berekend nadat de matrix getransponeerd wordt door t(ptm_matrix) 
Als Jaccard afstand = 0 dan komen PTM’s altijd samen voor. 
Als Jaccard afstand = 1 dan komen de PTM’s nooit samen voor. 
Op basis hiervan kan het samen voorkomen van PTM’s worden benaderd en gaan we deze afstanden reduceren tot 2 dimensies (k=2) zodat het daarna geclusterd kan worden met k-means.

#figure(
  image("ReportData/Picture17.png")
)

Deze coordinaten gaan we uitplotten als clusters om hun samen voorkomen te visualiseren.  We doen dit door k-means toe te passen: 

#figure(
  image("ReportData/Picture18.png")
)

We zorgen er eerst voor dat de random number generator van R een vaste seed heeft zodat de clusters reproduceerbaar zijn. Het aantal clusters kan zo nodig worden bijgesteld (k). Clusters worden gemaakt door het kmeans algoritme, deze minimaliseert de euclidische afstanden² in de gereduceerde dimensies. Met andere woorden gaan we simpelweg de dots die het dichst bij elkaar zijn groeperen. Het algormitme kiest ad random verschillende MDS_coord om te kijken wat de euclidische afstanden zijn met naburige coördinaten en zoekt dan de optimale clustering. 

#figure(
  image("ReportData/Picture19.png", width: 60%)
)

#figure(
  image("ReportData/Picture20.jpg"),
  caption: [Clusters zouden cellulaire/sub-cellulaire processen kunnen detecteren]

)

#figure(
  image("ReportData/Picture21.jpg"),
  caption: [Clusters zouden cellulaire/sub-cellulaire processen kunnen detecteren]
)
#v(10pt)
=== Statistical significance of PTM co-occurance
#v(10pt)
Om te weten of het samen voorkomen van PTM’s niet zomaar toevallig is kunnen we de Fisher exacte test gebruiken. Hiermee testen we of de associatie van 2 PTM’s doelmatig is en niet louter door kans, dit doen we met een p-waarde. We bouwen eerst een kruistabel op (contingency matrix) doormiddel van een geneste for-loop zodat we doorheen ptm_matrix kunnen itereren en de PTM’s definiëren in een nieuwe matrix. 

#figure(
  image("ReportData/Picture22.png")
)

We stellen dus een contingency matrix voor alle PTM paren, dankzij de binnenste loop met elementen j die begint vanaf i+1 zullen identieke PTM’s niet opgenomen worden in de kruistabel. De resulterende matrix heeft 2 rijen en 2 kolommen :

#figure(
  image("ReportData/table2.png")
)

Deze matrix laten we dan door de ingebouwde functie in R verwerken :

#figure(
  image("ReportData/Picture23.png", width: 80%)
)

fisher.test berekent de p-waarde als volgt : 

#figure(
  image("ReportData/picture24.png")
)

De nulhypothese van een Fisher exacte test is dat 2 PTM’s onafhankelijk zijn (of niet in interactie gaan). Vanaf welke p-waarde je iets significant vind is ieders persoonlijke keus. Onderstaand heb je 10 outputs als voorbeeld:

#figure(
  image("ReportData/Picture25.png")
)
De kans dat bijvoorbeeld carbamylatie en methylatie onafhankelijk voorkomen is bijvoorbeeld extreem klein. De kleinste p-waarden heb ik van boven laten verschijnen met order. 


#v(10pt)
== Opschaling in python
#v(10pt)
De PTM analyse werd ook geïmplementeerd in Python om het daar op te schalen voor alle .pepXML bestanden. Alle details van de code zijn terug te vinden in de Github repository, maar hier is een algemene overview van het proces:

1.	Opnieuw is er enkel interesse in betrouwbare target matches (expect <= 5%).
2.	Indien er een PTM-combinatie gematched kan worden met genoeg betrouwbaarheid, wordt dit opgeslagen in een lokale database.
3.	Wanneer de database gevuld is, kunnen de PTM matches per individu worden opgevraagd en geplot op een grafiek, waarna er eindelijk vergeleken kan worden tussen individuen. Alle plots staan eveneens op Github.

Het volgende kon uit de grafieken waargenomen worden:

1.	De gemeten frequenties per individu variëren, waarschijnlijk omdat niet elke betrouwbare target match betrouwbaar gematched kon worden met PTM’s.
2.	Carbamidomethyl, Carbamylation en Methylation zijn de meest voorkomende            PTM’s bij alle individuen over alle groepen.
3.	De volgende 6 PTM’s hadden de meest zichtbare variantie in de frequenties: Dimethylation, Dioxidation, Disulfideloss, Nitration, Nitrosylation and Oxidation.

Er werden geen significante verschillen waargenomen tussen de groepen.

#v(10pt)
== Substition finder
#v(10pt)
dit is een simpele tool als proof of concept om puntmutaties te vinden op basis van PSM’s. 
#v(10pt)
=== Core logic
#v(10pt)

uit de bemonsterde PSM’s wordt eerst de namen van het eiwit en de sequenties genomen  en verzameld in psm_data variabele. We kiezen hiervoor enkel de best gerangde hits van elk spectrum (de eerste). 

#figure(
  image("ReportData/Picture28.png")
)

Het is dan handiger om de variabele psm_data, die uit meerdere kleine dataframes bestaat, samen te voegen tot 1 grote dataframe, dit gaan we dan doen met rbind en PSM’s zonder hits gaan we uitfilteren : 

#figure(
  image("ReportData/Picture29.png")
)

#figure(
  image("ReportData/Picture30.png")
)

Nu hebben we een gebruiksvriendelijke variabele psm_df met alle eiwit ID’s en sequenties uit de sample chunk en kunnen we hieruit substituties zoeken ten opzichte van de FASTA. 

We gebruiken de grepl functie van R om pattern matches in strings te vinden om de, hiermee zoeken we naar het corresponderende eiwit ID van de FASTA die overeenkomt met de dataframe in psm_df. Als er meerdere matches zijn dan word er simpelweg de eerste genomen die gevonden word. De lengte (aantal letters) gaan we opslaan in de variabele peptide_len. 

Voor de effectieve vergelijking gaan we substrings van de FASTA sequenties uithalen die dezelfde lengte (window) hebben als de sequenties uit psm_df, deze zijn allemaal kandidaten om gealigneerd te worden. Deze substrings worden tegelijk dan ook een op een vergeleken met de sequentie van psm_df, de string word dus opgesplitst met de strsplit functie en alle letters (aminozuren) worden dan individueel vergeleken. Als er een verschil word gedetecteerd dan word deze als TRUE opgeslagen in de variabele diffs, dit aantal word bijgehouden en op een maximum van 1 gehouden (omdat de kans van 2 substituties op zulke korte sequenties klein zijn, als er toch 2 letters zouden verschillen zal het dan eerder een aparte sequentie zijn die toevallig gelijkend is). Dus als de peptide sequentie 1 aminozuur (1 letter) verschilt, dan word de plaats van die ongelijkheid aangeduidt en de verandering zelf ook opgeslagen (bv. K -> R ). 

voorbeeld outputs : 

#figure(
  image("ReportData/Picture31.png")
)

#v(10pt)
= PTM Analyze Conclusie
#v(10pt)
Dit onderzoek heeft zich gericht op de detectie en interpretatie van post-translationele modificaties (PTMs) in tumorstalen van diverse etnische groepen, met behulp van geavanceerde massaspectrometrie-gegevens en bio-informatica tools. Door middel van een gesubsamplepe pepXML-analyse in R werd een efficiënte methode ontwikkeld om PTM-patronen te identificeren en statistisch te valideren.

Enkele belangrijke bevindingen:
-  PTM-aanrijking en diagnostische waarde: De ontwikkelde PTM-matcher toonde aan dat bepaalde modificaties (zoals methylatie en carbamylatie) significant vaker samen voorkomen, wat mogelijk wijst op gedeelde regulatorische pathways in tumorbiologie.

-    Populatieverschillen: Hoewel verder onderzoek nodig is, suggereren de Jaccard-clustering en Fisher exacte tests dat sommige PTM-combinaties mogelijk vaker voorkomen in specifieke etnische groepen, wat zou kunnen wijzen op onderliggende genetische of omgevingsinvloeden.

-    Methodologische innovatie: Door subsampling en multidimensionale schaling (MDS) konden we rekenkracht beperken zonder significante informatieverlies, wat de deur opent voor schaalbare analyses van grote proteomics-datasets.

Deze inzichten vormen een stap richting gepersonaliseerde oncologie, waarbij PTM-profielen mogelijk kunnen bijdragen aan betere stratificatie van patiënten op basis van etnische achtergrond of tumorbiologie. Vervolgonderzoek zou zich moeten richten op functionele validatie van deze PTM-clusters en hun rol in ziektegerelateerde pathways.

Finaal gaan we kijken naar de PTM’s die specifiek het tumor suppressie pathway opmaken. Methylatie, phosphorylatie, sumoylatie, ubiquitinering en acetylering. Met behulp van een 2-way ANOVA en TukeyHSD post-hoc analyse kunnen we zien of binnen verschillende etnische groepen deze pathway verschilt in expressiegraad. Onderaan zijn de paarsgewijze vergelijkingen visueel duidelijk gemaakt met de boxplots. 

#figure(
  image("ReportData/Picture27.jpg")
)

de frequentie van de PTM’s verschilt niet significant tussen verschillende etnische groepen. Er is geen enkele p-waarde kleiner dan 0,05 gevonden. Het is duidelijk dat deze tumorsuppressie-pathway sterk geconserveerd is. Met andere woorden, deze moleculaire signaalcascade was al reeds aanwezig in een voorouderlijk organisme en werd dus overgeërfd naar de vroege Homo sapiens, dit verklaart de gelijkenis in PTM frequentie doorheen deze etnische groepen. 

Met dit project is niet alleen een reproduceerbare analytische pijplijn gecreëerd, maar ook een basis gelegd voor toekomstige studies naar proteomische diversiteit in verschillende populaties.

#pagebreak()


= Chan Ming Jan - Biologie
- Theorethische achtergrond rond proteometics en spectrometrie verduidelijkt.
- In R code PTM  subset analyzer geimplementeerd.
- In R code AA Sequence variant finder geimplementeerd.
- Conlusies trekken met ANOVA.
- Biologische steun voor de groep.
- Research resultaten tekst geleverd for finale report.

= Axel De Leeuw - Informatica
- Dataset door fragpipe laten proccesen.
- Research rond dataset, fragpipe en CalcUA (supercomputer) gebruik.
- Analyze in python rond Open vs Closed search geimplementeerd.
- Geresearched rond aminozuur substituties.
- PTM match graphs geleverd.
- Research resultaten geleverd for finale report.

= Anass Hamzaoui - Informatica
- Dataset door fragpipe laten proccesen.
- Research rond dataset, fragpipe en CalcUA (supercomputer) gebruik.
- Analyze in python rond Open vs Closed search geimplementeerd.
- Alle report material bij elkaar getrokken.
- Report in typst format in elkaar gestoken.
- Jan's R code geconvert naar python.

