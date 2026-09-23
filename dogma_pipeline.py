"""
dogma_pipeline.py — Pipeline del dogma central con Biopython.
Replicación (hebras complementarias) -> Transcripción (ARNm) -> Traducción (proteína).

Uso:
    python dogma_pipeline.py secuencia.fasta [--hebra codificante|molde] [--tabla 1] [--salida prefijo]
"""
import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

log = logging.getLogger("dogma")


def configurar_log(nivel=logging.INFO):
    if not log.handlers:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%H:%M:%S"))
        log.addHandler(h)
    log.setLevel(nivel)
    log.propagate = False


def corta(s, n=45):
    s = str(s)
    return s if len(s) <= n else f"{s[:n]}... ({len(s)} nt/aa)"


class PipelineDogmaCentral:
    def __init__(self, secuencia, nombre="secuencia", hebra="codificante", tabla=1):
        self.nombre, self.hebra, self.tabla = nombre, hebra, tabla
        self.adn = Seq(str(secuencia).upper().replace("U", "T"))
        self.resultados = {}
        invalidas = set(str(self.adn)) - set("ACGTN")
        if invalidas:
            raise ValueError(f"Caracteres no válidos en la secuencia: {invalidas}")
        log.info("=" * 70)
        log.info(f"Secuencia '{nombre}' cargada: {len(self.adn)} nt, hebra '{hebra}'")
        gc = 100 * (self.adn.count("G") + self.adn.count("C")) / max(len(self.adn), 1)
        log.info(f"Contenido GC: {gc:.1f} %")

    @classmethod
    def desde_fasta(cls, ruta, **kw):
        reg = next(SeqIO.parse(ruta, "fasta"))
        log.info(f"Leyendo FASTA '{ruta}': {reg.description}")
        return cls(reg.seq, nombre=reg.id, **kw)

    # ---------- PASO 1: REPLICACIÓN ----------
    def replicacion(self):
        log.info("-" * 70)
        log.info("PASO 1 · REPLICACIÓN (semiconservativa)")
        sup = self.adn if self.hebra == "codificante" else self.adn.reverse_complement()
        inf = sup.reverse_complement()  # hebra complementaria escrita 5'->3'
        log.info("Helicasa: separa la doble hélice en dos hebras molde")
        log.info("Primasa: coloca cebadores de ARN; ADN polimerasa: sintetiza 5'->3'")
        nueva_1 = sup.reverse_complement()   # sintetizada sobre la hebra superior
        nueva_2 = inf.reverse_complement()   # sintetizada sobre la hebra inferior
        log.info(f"Hebra parental 5'->3'           : {corta(sup)}")
        log.info(f"Hebra parental complementaria   : {corta(inf)}")
        log.info(f"Nueva hebra (molde = superior)  : {corta(nueva_1)}")
        log.info(f"Nueva hebra (molde = inferior)  : {corta(nueva_2)}")
        log.info("Ligasa: une los fragmentos de Okazaki de la hebra retrasada")
        ok = str(nueva_1) == str(inf) and str(nueva_2) == str(sup)
        log.info(f"Comprobación: moléculas hijas idénticas a la parental -> {ok}")
        self.codificante, self.molde = sup, inf
        self.resultados["replicacion"] = {"hija_1": (sup, nueva_1), "hija_2": (nueva_2, inf), "ok": ok}
        return self.resultados["replicacion"]

    # ---------- PASO 2: TRANSCRIPCIÓN ----------
    def transcripcion(self):
        if "replicacion" not in self.resultados:
            self.replicacion()
        log.info("-" * 70)
        log.info("PASO 2 · TRANSCRIPCIÓN")
        log.info("ARN polimerasa lee la cadena MOLDE 3'->5' y sintetiza ARN 5'->3' (U en lugar de T)")
        arnm = self.molde.reverse_complement().transcribe()
        log.info(f"Cadena molde 5'->3' : {corta(self.molde)}")
        log.info(f"ARNm 5'->3'         : {corta(arnm)}")
        log.info(f"Comprobación: ARNm == codificante con U -> {str(arnm) == str(self.codificante.transcribe())}")
        self.arnm = arnm
        self.resultados["transcripcion"] = arnm
        return arnm

    # ---------- búsqueda de ORF ----------
    def buscar_orf(self):
        tabla_stop = {"TAA", "TAG", "TGA", "UAA", "UAG", "UGA"}
        mejor = None
        s = str(self.arnm)
        for marco in range(3):
            i = marco
            while i + 3 <= len(s):
                if s[i:i+3] == "AUG":
                    j = i
                    while j + 3 <= len(s) and s[j:j+3] not in tabla_stop:
                        j += 3
                    fin = j + 3 if j + 3 <= len(s) else j
                    if mejor is None or fin - i > mejor[1] - mejor[0]:
                        mejor = (i, fin, marco, j + 3 <= len(s))
                    i = j
                i += 3
        if mejor is None:
            log.warning("No se ha encontrado ningún codón AUG: no hay ORF")
            return None
        ini, fin, marco, con_stop = mejor
        log.info(f"ORF más largo: posiciones {ini+1}-{fin}, marco {marco+1}, "
                 f"{'con' if con_stop else 'SIN'} codón de paro")
        return self.arnm[ini:fin]

    # ---------- PASO 3: TRADUCCIÓN ----------
    def traduccion(self):
        if "transcripcion" not in self.resultados:
            self.transcripcion()
        log.info("-" * 70)
        log.info("PASO 3 · TRADUCCIÓN")
        orf = self.buscar_orf()
        if orf is None:
            self.proteina = Seq("")
            return self.proteina
        log.info(f"Codón de inicio: {orf[:3]} (Met) -> el ribosoma ensambla y empieza en el sitio P")
        log.info(f"Codón de paro  : {orf[-3:]} -> factor de liberación termina la traducción")
        orf3 = orf[: len(orf) - len(orf) % 3]
        prot = orf3.translate(table=self.tabla, to_stop=True)
        log.info(f"Codones leídos : {len(orf3)//3}  ->  {len(prot)} aminoácidos")
        log.info(f"Proteína (N->C): {corta(prot, 60)}")
        if prot:
            pa = ProteinAnalysis(str(prot).replace("X", ""))
            log.info(f"Masa ≈ {pa.molecular_weight()/1000:.2f} kDa | pI ≈ {pa.isoelectric_point():.2f} | "
                     f"GRAVY = {pa.gravy():.3f}")
        self.proteina = prot
        self.resultados["traduccion"] = prot
        return prot

    def ejecutar(self):
        self.replicacion()
        self.transcripcion()
        self.traduccion()
        log.info("=" * 70)
        log.info("Pipeline completado")
        return self.resultados

    def guardar(self, prefijo):
        recs = [
            SeqRecord(self.codificante, id=f"{self.nombre}_codificante", description="ADN 5'->3'"),
            SeqRecord(self.molde, id=f"{self.nombre}_molde", description="ADN complementario 5'->3'"),
            SeqRecord(self.arnm, id=f"{self.nombre}_ARNm", description="ARNm 5'->3'"),
            SeqRecord(self.proteina, id=f"{self.nombre}_proteina", description="N->C"),
        ]
        ruta = f"{prefijo}_resultados.fasta"
        SeqIO.write(recs, ruta, "fasta")
        log.info(f"Resultados guardados en {ruta}")
        return ruta


def main():
    ap = argparse.ArgumentParser(description="Replicación -> Transcripción -> Traducción")
    ap.add_argument("fasta")
    ap.add_argument("--hebra", default="codificante", choices=["codificante", "molde"])
    ap.add_argument("--tabla", type=int, default=1, help="tabla de código genético NCBI")
    ap.add_argument("--salida", default=None, help="prefijo para guardar resultados FASTA")
    a = ap.parse_args()
    configurar_log()
    p = PipelineDogmaCentral.desde_fasta(a.fasta, hebra=a.hebra, tabla=a.tabla)
    p.ejecutar()
    if a.salida:
        p.guardar(a.salida)


if __name__ == "__main__":
    main()
