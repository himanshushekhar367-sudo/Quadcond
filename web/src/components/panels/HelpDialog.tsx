import { useRef } from 'react';
import { BookOpen, X } from 'lucide-react';
import { Button } from '@/components/ui/button';

export function HelpDialog() {
  const dialog = useRef<HTMLDialogElement>(null);
  return <>
    <Button variant="outline" onClick={() => dialog.current?.showModal()} aria-haspopup="dialog">
      <BookOpen className="h-4 w-4" aria-hidden /> Guide
    </Button>
    <dialog ref={dialog} aria-labelledby="guide-title" className="workbench-guide">
      <div className="flex items-center justify-between gap-4 border-b border-border pb-4">
        <h2 id="guide-title" className="text-xl font-semibold">Using QuadCond and AENNA-3D</h2>
        <Button variant="ghost" aria-label="Close guide" onClick={() => dialog.current?.close()}><X className="h-5 w-5" /></Button>
      </div>
      <div className="space-y-5 py-5 text-sm leading-relaxed text-fg-muted">
        <section><h3 className="mb-2 font-semibold text-fg">Start with a DNA sequence</h3>
          <p>Paste a sequence in the Sequence panel or choose an example from the library. Both strands are shown in the 5′ to 3′ direction. The entered strand is + and its reverse complement is −. G4 and i-motif models are routed to the appropriate strand, which is stated beside every result.</p>
        </section>
        <section><h3 className="mb-2 font-semibold text-fg">Read the evidence before the picture</h3>
          <p>Evidence cards report predicted melting temperature or transitional pH, their intervals and applicability. Select a card to highlight its schematic archetype. On mobile, select the card in Evidence, then open 3D. The picture is an illustration of a structural class, not an atomic structure prediction.</p>
          <p className="mt-2">Unsupported estimates are hidden until you request extrapolation. A refused prediction has no number and cannot be restored with this option. Genomic scores describe antibody peak associations, not folding probabilities or molecular occupancy.</p>
        </section>
        <section><h3 className="mb-2 font-semibold text-fg">Compare sequences and conditions</h3>
          <p>Open Compare to scan substitutions, sweep one condition or submit a batch. Mutation differences are exploratory; the spread across estimators is not a validated confidence interval for a mutation effect. Condition summaries use supported points only.</p>
          <p className="mt-2">Batch input accepts FASTA or one sequence per line. Empty records are reported as exclusions, and repeated names remain separate records. Download CSV or JSON to retain the input, conditions and model identity with the result.</p>
        </section>
        <section><h3 className="mb-2 font-semibold text-fg">Privacy and availability</h3>
          <p>The prediction service uses no accounts and does not persist submitted sequences. Requests are sent to the configured QuadCond service; deployment operators must describe any infrastructure logs and retention policy. A disconnected service provides schematic viewing with no predicted values.</p>
          <p className="mt-2">QuadCond backend: <a className="underline text-accent" href="/quadcond-license.txt" target="_blank" rel="noreferrer">MIT license</a>. <a className="underline text-accent" href="/quadcond-api-example.py" download>Download the Python API example</a>.</p>
        </section>
      </div>
    </dialog>
  </>;
}
