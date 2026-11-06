import type { ComponentProps } from "react";
import { cva, type VariantProps } from "class-variance-authority";
import { cn } from "@/lib/utils";

const badgeVariants = cva(
  "inline-flex items-center whitespace-nowrap rounded-[var(--radius-xs)] border px-1.5 py-1 text-[0.6875rem] font-medium leading-none tracking-wide",
  {
    variants: {
      variant: {
        default: "border-border bg-surface-2 text-fg-muted",
        accent: "border-transparent bg-accent text-accent-foreground",
        outline: "border-border text-fg-muted",
      },
    },
    defaultVariants: { variant: "default" },
  },
);

export function Badge({
  className,
  variant,
  ...props
}: ComponentProps<"span"> & VariantProps<typeof badgeVariants>) {
  return <span className={cn(badgeVariants({ variant }), className)} {...props} />;
}
