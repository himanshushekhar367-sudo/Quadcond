import type { EnsembleMember, Polymer, ProteinBinder, StructureKind } from "./types";

interface BinderRecord extends Omit<ProteinBinder, "score" | "rationale"> {
  motifs?: string[];
}

const CATALOG: BinderRecord[] = [
  {
    id: "hnrnpa1",
    name: "Heterogeneous nuclear ribonucleoprotein A1",
    gene: "HNRNPA1",
    uniprot: "P09651",
    kinds: ["g-quadruplex", "hairpin", "ss-coil"],
    polymers: ["DNA", "RNA"],
    organism: "Homo sapiens",
    motifs: ["GGG", "TTAGGG", "UUAGGG"],
    sequence:
      "MSKSESPKEPEQLRKLFIGGLSFETTDESLRSHFEQWGTLTDCVVMRDPNTKRSRGFGFVTYATVEEVDAAMNARPHKVDGRVVEPKRAVSREDSQRPGAHLTVKKIFVGGIKEDTEEHHLRDYFEQYGKIEVIEIMTDRGSGKKRGFAFVTFDDHDSVDKIVIQKYHTVNGHNCEVRKALSKQEMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQSSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF",
  },
  {
    id: "ncl",
    name: "Nucleolin",
    gene: "NCL",
    uniprot: "P19338",
    kinds: ["g-quadruplex", "hairpin"],
    polymers: ["DNA", "RNA"],
    organism: "Homo sapiens",
    motifs: ["GGGG", "MYC"],
    sequence:
      "MVKLAKAGKNQGDPKKMAPPPKEVEEDSEDEEMSEDEEDSSEDEEDSSEDEEEEEEEEEPVAKKGVKPKKATKKAEEAEEEPEEKKEKEMEKQKEKEPEKLKGKVEDTKEEKEEKKEEEKQEEEEDKEEEEEEEEKEEEDDEEDEDDEDEDEEDDDEDEDDEDDDEDDEEEEEEEEEEPVKEAPGKRKKEMAKQKAAPEAKKQKVEGTEPTTAFNLFVGNLNFNKSAPELKTGISDVFAKNDLAVVDVRTGTNRKFGYVDFESAEDLEKALELTGLKVFGNEIKLEKPKGKDSKKERDARTLLAKNLPYKVTQDELKEVFEDAAEIRLVSKDGKSKGIAYIEFKTEADAEKTFEEKQGTEIDGRSISLYYTGEKGQNQDYRGGKNSTWSGESKTLVLSNLSYSATEETLQEVFEKATFIKVPQNQNGKSKGYAFIEFASFEDAKAAVAKKGVEVDGKNVTLMFSAGESKTLVLGNLSYSATEETLQEVFEKATFIKVPQNQNGKSKGYAFIEFASFEDAKEALNSCNKREIEGRAIRLELQGPRGSPNARSQPSKTLFVKGLSEDTTEETLKESFDGSVRARIVTDRETGSSKGFGFVDFNSEEDAKAAKEAMEDGEIDGNKVTLDWAKPKGEGGFGGRGGGRGGFGGRGGGRGGRGGFGGRGRGGFGGRGGFRGGRGGGGDHKPQGKKTKFE",
  },
  {
    id: "dhx36",
    name: "ATP-dependent RNA helicase DHX36 (RHAU)",
    gene: "DHX36",
    uniprot: "Q9H2U1",
    kinds: ["g-quadruplex"],
    polymers: ["DNA", "RNA"],
    organism: "Homo sapiens",
    motifs: ["GGGT", "GGGU"],
    sequence:
      "MSNYTKLFDNLHQERLKDILKGEVVKFGGKGGKGGKGGKRGGRGRGRGAGNSETGRGGRGRGRGAGRGGGGGGSGGGGSGGGGKRGGRGQKGGSGKGGKGGKGGKGGKRGGRGGKGGKGGKGGSGKRGGRGGKGGSGKRGGSGKRGGRGQKGGSGKRGGKGGKGGSGKRGGKGGKRGGRGQKGESGKRGGRGQKGGSGKRGGRGQKGGSGKRGGKGGKRGGRGQKGESGKRGGRGQKGGSGKRGGRGQKGGSGKRGGKGGKRGGRGQKGESGKRGGRGQKGGSGKRGGRGQKGGSGKRGGKGGKRGGRGQKIIYTGAGTGKTTYTMMGILQQQLEKNPQIIICSPSRELANTQKLEVVTDGVALPPEKIGVILNEMKFNPDLIRGKVIVFTQTKKEADELTRLGCHVVVATPGRLLDHLENTVGFSLIQKDLLKLRKSLKPGGKPDPKLVKESHPDVVVYELIPEYISKLLLNPY",
  },
  {
    id: "cnbp",
    name: "Cellular nucleic acid-binding protein",
    gene: "CNBP",
    uniprot: "P62633",
    kinds: ["g-quadruplex", "ss-coil"],
    polymers: ["DNA", "RNA"],
    organism: "Homo sapiens",
    sequence:
      "MSSNECFKCGRSGHWARECPTGGGRGRGMRSRGRGFQKRQKVGHPPCPVERRCGECGGSGHWAEACNTARGNGRRFHRRKKQCRYSCPTEGRCHGCGRGGHWAAECTAQARGRGHAAEHPRPDDQECCLICNKKAGPWEKTCPETPKRVRGGRGRGRRGAGNSET",
  },
  {
    id: "fmr1",
    name: "Fragile X messenger ribonucleoprotein 1",
    gene: "FMR1",
    uniprot: "Q06787",
    kinds: ["g-quadruplex", "hairpin"],
    polymers: ["RNA"],
    organism: "Homo sapiens",
    motifs: ["GGGU", "G-quadruplex"],
    sequence:
      "MEELVVEVRGSNGAFYKAFKNDSEMTVRAFEMEAAIRQEELAPVNDFNHNYKFAAESIKQPRNLNKKKSEADLDKLNEDQKAGEQAVHQEFVHELQLRANNSRIKNNFFNDKTSRHHHLAAEQLYNSVVQASTKVQQYKNVQKEADRCLKDLPEDTIPQFTDLKLSVLENYQSVHSQTAKQAAAYSDHNDLSDNRQNTKAFSNVTQAANEKEEFVHLERNLQKEQNDLQELNQSNKADLEELQKLEQERDQYAEKLQDFERQKQKLEELNKSNEELKLEKHSN",
  },
  {
    id: "pot1",
    name: "Protection of telomeres protein 1",
    gene: "POT1",
    uniprot: "Q9NUX5",
    kinds: ["g-quadruplex", "ss-coil"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    motifs: ["TTAGGG", "GGGTTA"],
    sequence:
      "MSLVPATNYIYTPLNQLKGGTIVNVYGVVKFFKPPYLSKGTDYCSVVTIVDQTNVKLTCLLFSGNYEALPIIYKNGDIVRFHRLSKVHYTTVNTGKRQKGAEKTALSTRFSIKNRYYVPLKDIWMDICNVTLTEAKFALQSTGTKIVYLHCTKCNLKDTNLIQKHAQIHKHQGIYVVDQRKSTNTTNIQKLTLYECFQGQINLNSLIAELQKKQHQSYLNNKK",
  },
  {
    id: "blm",
    name: "Bloom syndrome protein",
    gene: "BLM",
    uniprot: "P54132",
    kinds: ["g-quadruplex", "cruciform", "duplex"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    sequence:
      "MAAVLEANAEEPKCWNPQNAKPFYFNYAAIQDLPIEKLHVDFVCADNSVILTNYLIIDWSWNKVSQTKQVQLLQKVDPDYYLNKLKAEKQVSLRSLFNKYTDKIIQDGITKLKDLYTDKENAIPQVVGKVQNKVLKGPKAAKCVNLGKEILKMSNVCQSNLEAINAALNKDQVAITKLQKELTQLKNDCSNFNSQLTSYQKLSKEDLQNLKNKDFGCLKKLLSKFEGLTDSKSYQELQKKLEALKNKDF",
  },
  {
    id: "hnrnpk",
    name: "Heterogeneous nuclear ribonucleoprotein K",
    gene: "HNRNPK",
    uniprot: "P61978",
    kinds: ["i-motif", "ss-coil", "hairpin"],
    polymers: ["DNA", "RNA"],
    organism: "Homo sapiens",
    motifs: ["CCC", "CCT"],
    sequence:
      "METEQPEETFPNTETNGEFGKRPAEDMEEEQAFKRSRNTDEMVELRILLQSKNAGAVIGKGGKNIKALRTDYNASVSVPDSSGPERILSISADIETIGEILKKIIPTLEGLLETKFSPKGQRGLSFEAAGKINTLIINGRPIRNITELNKPPDTPGRRPITITGTQDQIQNAQYLLQNSVKQYSGKFFCGRKEPDPEKGKRKRRESEKEDDEEPPLSPVLCTSLPPSTTDAHSPEYSNAGAVKPKAQPADAGSKPLNHGDDDGFGDSDEN",
  },
  {
    id: "bcl6",
    name: "B-cell lymphoma 6 protein",
    gene: "BCL6",
    uniprot: "P41182",
    kinds: ["i-motif", "duplex"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    motifs: ["CCC"],
    sequence:
      "MASPADSCIQFTRHASDVLLNLNRLRSRDILTDVVIVVSREQFRAHKTVLMACSGLFYSIFTDQLKCNLSVINLDPEINPEGFCILLDFMYTSRLNLREGNIMAVMATAMYLQMEHVVDTCRKFIKASEAEMKSEMILKHKKKCFQTFNGTQALHIGCRNSIRAPASPDSSSSSSSSSSSSSSSSSSSSSPSPPPLPIE",
  },
  {
    id: "lin28a",
    name: "Protein lin-28 homolog A",
    gene: "LIN28A",
    uniprot: "Q9H9Z2",
    kinds: ["hairpin", "ss-coil"],
    polymers: ["RNA"],
    organism: "Homo sapiens",
    sequence:
      "MGSVSNQQFAGGCAKAAEEAPEEAPEDAARAADEPQLLHGAGICKWFNVRMGFGFLSMTARAGVALDPPVDVFVHQSKLHMEGFRSLKEGEAVEFTFKKSAKGLESIRVTGPGGVFCIGSERRPKGKSMQKRRSKGDRCYNCGGLDHHAKECKLPPQPKKCHFCQSISHMVASCPLKAQQGPSAQGKPTYFREEEEEIHFSPT",
  },
  {
    id: "drosha",
    name: "Ribonuclease 3 (DROSHA)",
    gene: "DROSHA",
    uniprot: "Q9NRR4",
    kinds: ["hairpin"],
    polymers: ["RNA"],
    organism: "Homo sapiens",
    sequence:
      "MPAGVRLNRAREAKGRGSRGSGGGGGGGRGRGRGQGQGQGQGQGQGQGQNPSRGGGGGGGGGGGGGGGGGGSGGGGSGGGGSGGGGSGGGGSGGGGSKLNEKAAGLGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGSGGGGSGGGGSGGGGSGGGGSGGGGSKPDPELQKLLEQLEKQLEKAEKQLKQQLEKAEKQLKQQLEKA",
  },
  {
    id: "ssb",
    name: "Single-stranded DNA-binding protein",
    gene: "SSBP1",
    uniprot: "Q04837",
    kinds: ["ss-coil", "hairpin"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    sequence:
      "MFRRPVLQVLRQFVRHEMLRKMLLNHYKVWIKKIGEAKKFGVEAAKTRDTLEEVQKMNLEKISSLRNLSVAAKKTKAQSLAKTAARKQRPETVAATQTTAQQKRPETVAATQTTAQQKR",
  },
  {
    id: "tp53",
    name: "Cellular tumor antigen p53",
    gene: "TP53",
    uniprot: "P04637",
    kinds: ["duplex", "cruciform"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    sequence:
      "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD",
  },
  {
    id: "ctcf",
    name: "Transcriptional repressor CTCF",
    gene: "CTCF",
    uniprot: "P49711",
    kinds: ["duplex", "cruciform"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    sequence:
      "MEGDAVEAIVEESETFIKGKERKTYQCKQCSRTFTRRHHLLRHNRIHTGEKPYQCEYCSKTFRTRRHHLLRHNRIHTGEKPYQCDYCSKTFRTRRHHLLRHNRIHTGEKPYQCDYCSKTFSTRKHHLLRHNRIHTGEKPYSCDYCSKTFSTRKHHLLRHNRIHTGEKPYSCDYCSKTFSTRKHHLLRHNRIHTGEKPYSCDYCSKTFSTRKHHLLRHNRIHTGEKPYSCNYCSKTFSTRKHHLLRHNRIH",
  },
  {
    id: "recq1",
    name: "ATP-dependent DNA helicase Q1",
    gene: "RECQL",
    uniprot: "P46063",
    kinds: ["cruciform", "duplex", "hairpin"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    sequence:
      "MASVSALSTFDLTEAWQEEVQELQKKLSLEQKELQQLQEELEKQQKELEKQQKELEKQQKELEKLQKELEKLQKELEKLQKELEKLQKELEKLQKELEKLQKELEKLQKELEKLQKELEKLQKELTQLKNDCSNFNSQLTSYQKLSKEDLQNLKNKDFGCLK",
  },
  {
    id: "hmgb1",
    name: "High mobility group protein B1",
    gene: "HMGB1",
    uniprot: "P09429",
    kinds: ["cruciform", "duplex", "triplex"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    sequence:
      "MGKGDPKKPRGKMSSYAFFVQTCREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGETKKKFKDPNAPKRPPSAFFLFCSEYRPKIKGEHPGLSIGDVAKKLGEMWNNTAADDKQPYEKKAAKLKEKYEKDIAAYRAKGKPDAAKKGVVKAEKSKKKKEEEEDEEDEEDEEEEEDEEDEDEEEDDDDE",
  },
  {
    id: "pc4",
    name: "Activated RNA polymerase II transcriptional coactivator p15",
    gene: "SUB1",
    uniprot: "P53999",
    kinds: ["ss-coil", "triplex"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    sequence:
      "MPKSKELVSSSSSGSDSDSEVDKKLKRKKQVAPEKPVKKQKTGETSRALSSSKQSSSSDDRVKKEKKVKKPEEEEVAQKLAEKKEEEMKKEEEMKKAEKEEAEKAEKEEAEKAEKEEAEK",
  },
  {
    id: "asf1a",
    name: "Histone chaperone ASF1A",
    gene: "ASF1A",
    uniprot: "Q9Y294",
    kinds: ["duplex"],
    polymers: ["DNA"],
    organism: "Homo sapiens",
    sequence:
      "MAKVSVLNVAVLENPSPFHVSQFLFESHLEKFSAELEYFHKLYLDNLELTIPNAVAKDIINYVQDVKTGKCAKLEKNVQENTLKEVNKILELFPQELVNNFVQEYKDELLQKLEKNVQENTLKEVNKILELFPQELVNNFVQEYKDELLQK",
  },
];

function motifHit(seq: string, motifs?: string[]): number {
  if (!motifs?.length) return 0;
  let n = 0;
  for (const m of motifs) {
    if (seq.includes(m)) n += 1;
  }
  return n;
}

export function bindersFor(
  seq: string,
  polymer: Polymer,
  ensemble: EnsembleMember[],
): ProteinBinder[] {
  const top = ensemble.filter((m) => m.probability > 0.04);
  const kinds = new Set(top.map((m) => m.kind));
  const out: ProteinBinder[] = [];

  for (const rec of CATALOG) {
    if (!rec.polymers.includes(polymer)) continue;
    const overlap = rec.kinds.filter((k) => kinds.has(k));
    if (!overlap.length) continue;
    const occ = overlap.reduce((a, k) => {
      const m = ensemble.find((e) => e.kind === k);
      return a + (m?.probability ?? 0);
    }, 0);
    const motif = motifHit(seq, rec.motifs);
    const score = Math.min(0.99, 0.35 + 0.5 * occ + 0.08 * motif);
    const rationale = `Prefers ${overlap.map((k) => k).join(", ")} on ${polymer}. Occupancy-weighted match ${Math.round(occ * 100)}%${motif ? `; sequence motifs (${rec.motifs?.join(", ")}) present.` : "."}`;
    out.push({ ...rec, score: Math.round(score * 100) / 100, rationale });
  }

  return out.sort((a, b) => b.score - a.score);
}

export function kindHasBinders(kind: StructureKind): boolean {
  return CATALOG.some((c) => c.kinds.includes(kind));
}
