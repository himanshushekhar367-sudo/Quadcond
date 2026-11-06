/**
 * Text in a number box, read as a value only when it is one the field accepts.
 *
 * Kept out of the component so the rule is testable: a draft such as "2" on
 * the way to "25" in a field whose minimum is 4 is not a value yet, and must
 * neither be committed nor thrown away.
 */
export function parseDraft(text: string, min: number, max: number): number | null {
  const t = text.trim();
  if (!t) return null;
  const v = Number(t);
  return Number.isFinite(v) && v >= min && v <= max ? v : null;
}
