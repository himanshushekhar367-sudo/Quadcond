import { useEffect, useMemo, useRef } from "react";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import * as THREE from "three";
import type { Nucleotide } from "@/lib/na/types";
import type { GeometryOnlyMember } from "@/lib/na/geometry-only";

import { useNA } from "@/lib/na/store";

const BASE_COLOR: Record<string, number> = {
  A: 0x8aa4c8,
  T: 0xc4b49a,
  U: 0xc4b49a,
  G: 0x88b49a,
  C: 0xc49292,
  N: 0x9aa0a6,
};

/**
 * A line is built as a THREE object and mounted through `<primitive>` rather
 * than written as `<line>`.
 *
 * React 19's JSX types resolve the intrinsic `line` to the *SVG* element, so
 * `<line geometry={...}>` typechecks against SVGLineElement and fails -- the
 * three.js element of the same name loses the name collision. `<primitive>`
 * sidesteps it, and disposal is explicit because the geometry and material are
 * now ours to own.
 */
function Line3D({
  points,
  color,
  opacity,
}: {
  points: THREE.Vector3[];
  color: number;
  opacity: number;
}) {
  const line = useMemo(() => {
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({ color, transparent: true, opacity });
    return new THREE.Line(geometry, material);
  }, [points, color, opacity]);

  useEffect(
    () => () => {
      line.geometry.dispose();
      (line.material as THREE.Material).dispose();
    },
    [line],
  );

  return <primitive object={line} />;
}

function StrandLine({ points, licorice }: { points: THREE.Vector3[]; licorice: boolean }) {
  return <Line3D points={points} color={0xd8dce0} opacity={licorice ? 0.8 : 0.55} />;
}

function TetradBondLine({ from, to, s }: { from: Nucleotide; to: Nucleotide; s: number }) {
  const points = useMemo(() => [
    new THREE.Vector3(from.x * s, from.y * s, from.z * s),
    new THREE.Vector3(to.x * s, to.y * s, to.z * s),
  ], [from, to, s]);
  return <Line3D points={points} color={0x22d3ee} opacity={0.6} />;
}

function StructureGroup({ member, offsetX }: { member: GeometryOnlyMember; offsetX: number }) {
  const visMode = useNA((s) => s.visMode);
  const s = 0.18;
  const strands = useMemo(() => {
    return member.strands
      .map(ids => {
        const pts = ids
          .map(i => member.nucleotides[i])
          .filter((n): n is Nucleotide => Boolean(n))
          .map(n => new THREE.Vector3(n.x * s, n.y * s, n.z * s));
        return pts.length >= 2 ? pts : null;
      })
      .filter((pts): pts is THREE.Vector3[] => pts !== null);
  }, [member]);

  const ionY = useMemo(() => {
    if (member.kind === "g-quadruplex") {
      const ys = member.nucleotides.filter(n => n.role === "tetrad").map(n => n.y * s);
      if (ys.length) return (Math.min(...ys) + Math.max(...ys)) / 2;
    }
    return null;
  }, [member]);

  return (
    <group position={[offsetX, 0, 0]}>
      {member.nucleotides.map((n, i) => {
        const hex = BASE_COLOR[n.base] ?? BASE_COLOR.N!;
        const r = n.paired ? 0.42 : 0.34;
        const scale = visMode === 'space-filling' ? 0.9 : (visMode === 'licorice' ? 0.12 : r);
        return (
          <mesh key={i} position={[n.x * s, n.y * s, n.z * s]} scale={scale}>
            <sphereGeometry args={[1, 18, 18]} />
            <meshStandardMaterial color={hex} roughness={visMode === 'space-filling' ? 0.8 : 0.38} metalness={0.08} emissive={hex} emissiveIntensity={0.12} />
          </mesh>
        );
      })}
      {visMode !== 'space-filling' && strands.map((pts, i) => (
        <StrandLine key={i} points={pts} licorice={visMode === 'licorice'} />
      ))}
      {/* Tetrad Hoogsteen hydrogen bonds (cyan lines connecting 4 corners) */}
      {visMode !== 'space-filling' && member.tetradBonds?.map(([a, b], i) => {
        const ntA = member.nucleotides[a];
        const ntB = member.nucleotides[b];
        if (!ntA || !ntB) return null;
        return <TetradBondLine key={`tb-${i}`} from={ntA} to={ntB} s={s} />;
      })}
      {/* Central coordinating ion (K+/Na+) */}
      {ionY !== null && (
        <mesh position={[0, ionY, 0]} scale={visMode === 'space-filling' ? 0.35 : 0.15}>
          <sphereGeometry args={[1, 18, 18]} />
          <meshStandardMaterial color={0x9333ea} roughness={0.38} emissive={0x9333ea} emissiveIntensity={0.3} />
        </mesh>
      )}
    </group>
  );
}

function SceneContents({ members, overlay, autoRotate }: { members: GeometryOnlyMember[]; overlay: boolean; autoRotate: boolean }) {
  const rootRef = useRef<THREE.Group>(null);
  
  useFrame((_, delta) => {
    if (autoRotate && rootRef.current) {
      rootRef.current.rotation.y += 0.18 * delta;
    }
  });

  // Overlay draws the archetypes it was handed, in the order it was handed
  // them. It used to filter on `m.probability > 0.08` -- which ranked the
  // pictures by a locally computed Boltzmann weight, exactly the thing the
  // viewer must not do. When that field was removed the comparison became
  // `undefined > 0.08`, silently false, and the overlay drew nothing at all.
  const shown = overlay ? members.slice(0, 3) : members.slice(0, 1);

  return (
    <group ref={rootRef}>
      {shown.map((m, i) => (
        <StructureGroup key={m.id || i} member={m} offsetX={overlay ? (i - 1) * 7.5 : 0} />
      ))}
    </group>
  );
}

export function NaScene({
  members,
  overlay,
  autoRotate,
}: {
  members: GeometryOnlyMember[];
  overlay: boolean;
  autoRotate: boolean;
}) {
  // Overlay draws the archetypes it was handed, in the order it was handed
  // them. It used to filter on `m.probability > 0.08` -- which ranked the
  // pictures by a locally computed Boltzmann weight, exactly the thing the
  // viewer must not do. When that field was removed the comparison became
  // `undefined > 0.08`, silently false, and the overlay drew nothing at all.
  const shown = overlay ? members.slice(0, 3) : members.slice(0, 1);

  if (!shown.length) {
    return (
      <div className="flex h-full items-center justify-center text-sm text-fg-muted">No structure</div>
    );
  }

  return (
    <div className="h-full w-full" style={{ touchAction: "none" }}>
      <Canvas camera={{ position: [12, 8, 16], fov: 42 }} gl={{ antialias: true }}>
        <color attach="background" args={[0xf4f8f6]} />
        <hemisphereLight args={[0xdfe4e8, 0x1a1e22, 0.7]} />
        <directionalLight position={[8, 12, 6]} intensity={1.15} color={0xf2f4f6} />
        <ambientLight intensity={0.18} color={0xffffff} />
        
        <SceneContents members={members} overlay={overlay} autoRotate={autoRotate} />
        
        <OrbitControls
          enableDamping
          dampingFactor={0.08}
          minDistance={6}
          maxDistance={48}
        />
      </Canvas>
    </div>
  );
}
