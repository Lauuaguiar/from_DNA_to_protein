# Del ADN a la proteína: replicación, transcripción y traducción con Biopython

**Autoras:** Amai Suárez Navarro y Laura Aguia Pérez
**Asignatura:** Bioinformática
**Fecha:** septiembre de 2026
**Repositorio:** [ENLACE AL REPOSITORIO]

---

## Resumen

En este trabajo resolvemos seis ejercicios sobre el dogma central de la biología molecular. Primero los hicimos a mano y después los programamos en Python con Biopython para comprobar que coincidían. Los cuatro primeros van por partes (replicación, transcripción, traducción y splicing alternativo). El quinto trata la relación entre secuencia y estructura de una proteína, y el sexto junta todo en un pipeline que parte de un FASTA real, el de la insulina humana, y llega hasta la proteína. Todo el código está en un notebook de Jupyter y en el script `dogma_pipeline.py`. En todos los casos el resultado manual y el automático fueron iguales, lo cual no es tan obvio como parece: los fallos que tuvimos por el camino casi siempre venían de confundir la orientación de las hebras.

## 1. Introducción

El dogma central describe cómo fluye la información genética: el ADN se copia a sí mismo, se transcribe a ARN mensajero y ese ARNm se traduce a proteína. Sobre el papel son tres pasos sencillos. A la hora de programarlos aparecen detalles que en clase se pasan por alto, como que un FASTA siempre se escribe 5'→3' o que la cadena que "se lee" no es la que tiene la misma secuencia que el ARNm.

El objetivo era entender cada etapa haciendo los cálculos a mano y luego ver si el ordenador daba lo mismo. También queríamos pensar qué pasa cuando algo sale mal (una base mal copiada, un codón de inicio mutado, un exón que se salta) porque es ahí donde la teoría se vuelve interesante.

## 2. Material y métodos

Usamos Python 3 con **Biopython 1.88** (Cock et al., 2009) y pandas. Las funciones principales fueron `Seq.complement()`, `reverse_complement()`, `transcribe()` y `translate()`, este último con las opciones `to_stop` y `cds`. Para las propiedades de las proteínas usamos `ProteinAnalysis` y la escala de hidropatía de Kyte y Doolittle (1982).

Las secuencias salieron de tres bases de datos públicas:

- NCBI Nucleotide, a través de `Bio.Entrez`: CDS de la insulina humana (*INS*, RefSeq NM_000207).
- Ensembl, mediante su API REST: transcritos del gen *FGFR2*.
- RCSB PDB: estructura de la ubiquitina humana (1UBQ).

El pipeline del ejercicio 6 es una clase de Python que ejecuta los tres pasos en orden y va escribiendo con `logging` lo que hace en cada momento. Se puede usar desde el notebook o desde la terminal.

## 3. Resultados y discusión

### 3.1 Replicación del ADN

Partimos de la doble hebra `5'-ATG CCG TTA GCT-3'` / `3'-TAC GGC AAT CGA-5'`. Al separarse, cada hebra sirve de molde para una nueva complementaria y antiparalela. Salen dos moléculas hijas idénticas a la original, cada una con una hebra vieja y una nueva. Esto es la replicación semiconservativa.

La hebra nueva que se forma sobre la superior es `3'-TAC GGC AAT CGA-5'`, que escrita en el sentido habitual queda `5'-AGC TAA CGG CAT-3'`. Biopython dio exactamente lo mismo con `complement()` y `reverse_complement()`.

En el proceso participan varias enzimas. La helicasa abre la doble hélice rompiendo los puentes de hidrógeno. La primasa pone un cebador corto de ARN, porque la ADN polimerasa no sabe empezar desde cero y necesita un extremo 3'-OH. La ADN polimerasa alarga siempre en sentido 5'→3': en la hebra líder lo hace de forma continua y en la retrasada a trozos (fragmentos de Okazaki). Al final, la ligasa une esos trozos. En el notebook hicimos una simulación de la horquilla con fragmentos de 4 nucleótidos; es un modelo de juguete, pero ayuda bastante a ver por qué una de las hebras tiene que ir "a saltos".

Si la polimerasa mete una base equivocada y ni su actividad correctora ni la reparación de desapareamientos lo arreglan, el error se fija en la siguiente ronda. Lo simulamos: pusimos una T donde tocaba una G y, tras una segunda replicación, las dos hebras llevaban ya el cambio. El codón CCG pasó a CAG y la prolina se convirtió en glutamina (`MPLA` → `MQLA`). A partir de ahí es una mutación heredable.

### 3.2 Transcripción

En la secuencia `5'-ATG CCT GAA TGC-3'` / `3'-TAC GGA CTT ACG-5'`, la cadena molde es la inferior. La ARN polimerasa la lee de 3' a 5' y fabrica el ARNm de 5' a 3': `5'-AUG CCU GAA UGC-3'`. La superior es la codificante, y se reconoce porque empieza por ATG y coincide con el ARNm salvo por las U.

El promotor no aparece en el fragmento. Estaría antes del inicio de la transcripción, hacia el extremo 5' de la codificante, y es donde se une la ARN polimerasa (en eucariotas con la caja TATA hacia la posición −25). La región codificante empieza en el ATG, y en este caso son las 12 bases que tenemos.

El experimento de cambiar la orientación fue lo más útil del ejercicio. Cuando el script sabe qué hebra es y la secuencia está en 5'→3', sale Met-Pro-Glu-Cys. Si le damos la molde diciéndole que es la codificante sale `GCA UUC AGG CAU` (Ala-Phe-Arg-His), y si pegamos la molde tal y como está escrita en el enunciado (3'→5') sale `CGU AAG UCC GUA`. Ninguna de las dos tiene sentido biológico. El ordenador no avisa, simplemente hace lo que le pides.

### 3.3 Traducción y mutaciones

El transcrito `5'-AUG UAU GCU UAA-3'` tiene el codón de inicio AUG y el de paro UAA. Se traduce como Met-Tyr-Ala, con Met en el extremo N. `translate(cds=True)` lo valida y devuelve `MYA`.

La mutación AUG → GUG dio más juego del esperado. Con la tabla estándar, GUG se lee como valina y Biopython rechaza la secuencia como CDS. Con la tabla bacteriana (tabla 11) sí la acepta, porque en procariotas GUG puede iniciar con metionina. En eucariotas el ribosoma buscaría el siguiente AUG, y resulta que hay uno escondido en `GUG U|AUG| CU UAA`, en otro marco de lectura. Desde ahí saldría Met-Leu y ya no hay codón de paro en fase. No lo habíamos visto haciéndolo a mano.

Si desaparece el codón de paro (UAA → CAA), el ribosoma sigue leyendo la 3' UTR. Con una UTR inventada pasamos de `MYA` a `MYAQAGF`. Si tampoco hay paro ahí, llega a la cola poli-A y va añadiendo lisinas, y en la célula ese ARNm acabaría degradándose (*non-stop decay*).

### 3.4 Splicing alternativo

Con un gen de cinco exones diseñamos cuatro isoformas: la completa (1-2-3-4-5), 1-2-4-5, 1-3-5 y 1-2-5. Para que se viera el problema del marco de lectura, los exones 3 y 4 tienen 11 y 13 nucleótidos, que no son múltiplos de 3 pero suman 24.

| Isoforma | Exones | Nt omitidos | ¿Marco conservado? | Proteína (aa) |
|---|---|---|---|---|
| A | 1-2-3-4-5 | 0 | Sí | 20 |
| B | 1-2-4-5 | 11 | No | 18 |
| C | 1-3-5 | 25 | No | 15 |
| D | 1-2-5 | 24 | Sí | 12 |

La isoforma D pierde 8 aminoácidos pero el resto de la proteína se mantiene igual. En B y C todo lo que va después del salto cambia. Es decir, lo que importa no es cuánto se quita sino si el número de nucleótidos es múltiplo de 3.

Esto permite sacar varias proteínas de un solo gen. Con unos 20 000 genes codificantes, el ser humano tiene un proteoma bastante mayor, y más del 90 % de los genes con varios exones sufren splicing alternativo. En Ensembl consultamos *FGFR2*, que tiene decenas de transcritos. El caso más conocido son las isoformas IIIb y IIIc, que se diferencian en un exón mutuamente excluyente del tercer dominio de tipo inmunoglobulina. Eso cambia el sitio de unión: la IIIb, epitelial, reconoce FGF7 y FGF10, y la IIIc, mesenquimal, otros FGF distintos. Un cambio de unos cincuenta aminoácidos decide con qué ligando trabaja el receptor.

### 3.5 Proteínas: secuencia, estructura y función

En el péptido Met-Ile-Ser-Gly-Val-Lys-His el extremo N es la Met (grupo amino libre) y el C la His (grupo carboxilo libre).

El orden de los aminoácidos determina cómo se pliega la proteína (Anfinsen, 1973). Los residuos hidrofóbicos tienden a esconderse en el interior y los polares quedan hacia el agua. Si una mutación cambia un hidrofóbico interno por uno hidrofílico, se mete una carga o un grupo polar donde no cabe y el núcleo se desestabiliza. Lo comprobamos con el GRAVY: el péptido original da 0,329 y al cambiar Val5 por Lys baja a −0,829, con el pI pasando de 8,5 a 10. En una proteína real eso podría romper el plegamiento o hacer que se agregue. El caso contrario también existe: en la anemia falciforme un Glu de la superficie de la hemoglobina pasa a Val, que es hidrofóbica, y las moléculas se pegan entre sí.

Para ver estructuras reales descargamos la ubiquitina (1UBQ; Vijay-Kumar et al., 1987) y usamos los registros HELIX y SHEET del fichero PDB para marcar cada residuo. Tiene una hélice α y una lámina β en solo 76 aminoácidos. Nos pareció un buen ejemplo porque es pequeña y se ve todo de un vistazo. Una sola prolina en mitad de la hélice bastaría para romperla, ya que no puede formar el puente de hidrógeno del esqueleto.

### 3.6 Pipeline integrador

El pipeline lee la CDS de la insulina humana (NM_000207, 333 nt), genera las dos hebras complementarias y comprueba que las hijas son idénticas. Después transcribe la molde, busca el ORF más largo y lo traduce. Obtiene la preproinsulina de 110 aminoácidos, que empieza por el péptido señal MALWMRLL... En cada paso escribe un mensaje con lo que hace, las secuencias y alguna comprobación, y al final guarda todo en un FASTA.

Para contestar qué punto es más vulnerable, metimos 2000 mutaciones puntuales aleatorias en la CDS. Salieron un [__] % silenciosas, un [__] % de sentido erróneo y un [__] % que creaban un paro prematuro. Estos porcentajes dependen del código genético y no de la etapa en la que ocurra el error. Lo que sí depende de la etapa es cuánto dura el daño:

- **En la replicación**, el error queda en el ADN, pasa a todos los ARNm y a todas las proteínas de la célula, y además se hereda.
- **En la transcripción**, afecta solo a un ARNm. En nuestra simulación, 1 de cada 10 proteínas salía mal frente a 10 de 10 con el error en el ADN.
- **En la traducción**, se estropea una única molécula de proteína.

Nuestra conclusión es que la replicación es el punto más delicado. Lo curioso es que la ADN polimerasa es la más precisa de las tres (en torno a un error por cada 10⁹-10¹⁰ bases contando la reparación, frente a unos 10⁻⁴ del ribosoma). Pero sus errores son los únicos que se quedan. Dentro de la secuencia, lo más peligroso son los cambios en el codón de inicio, los que generan paros prematuros, las inserciones o deleciones que desplazan el marco y las mutaciones en las señales de splicing.

## 4. Conclusiones

Programar el dogma central obliga a ser más preciso que hacerlo a mano. La orientación de las hebras fue la principal fuente de errores, y el ordenador no los detecta por sí solo. Los ejercicios de mutaciones fueron los que más nos enseñaron, porque aparecieron cosas que no esperábamos (el AUG escondido en otro marco, por ejemplo). El pipeline final funciona con cualquier FASTA y se puede reutilizar para otras secuencias.

## Referencias

- Alberts, B. et al. (2022). *Molecular Biology of the Cell* (7.ª ed.). W. W. Norton.
- Anfinsen, C. B. (1973). Principles that govern the folding of protein chains. *Science*, 181, 223-230.
- Cock, P. J. A. et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.
- Kyte, J. y Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology*, 157, 105-132.
- Vijay-Kumar, S., Bugg, C. E. y Cook, W. J. (1987). Structure of ubiquitin refined at 1.8 Å resolution. *Journal of Molecular Biology*, 194, 531-544.
- Ensembl: https://www.ensembl.org · NCBI Nucleotide: https://www.ncbi.nlm.nih.gov/nuccore · RCSB PDB: https://www.rcsb.org

---

## Cómo usar el repositorio

```
pip install biopython pandas jupyter
jupyter notebook bioinformatica_dogma_central.ipynb
```

El pipeline también se puede lanzar desde la terminal:

```
python dogma_pipeline.py INS_cds.fasta --salida INS
```

Antes de descargar de NCBI hay que poner un correo en `Entrez.email` (en el ejercicio 6 del notebook).

Contenido:

- `bioinformatica_dogma_central.ipynb`: notebook con los seis ejercicios.
- `dogma_pipeline.py`: pipeline replicación → transcripción → traducción.
- `informe.md`: este mismo informe.
- `*.fasta`: secuencias de entrada y resultados.
