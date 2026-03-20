/* tslint:disable */
/* eslint-disable */

/**
 * Find commensurate twist angles for a bilayer system.
 *
 * # Arguments
 * * `radius`           – truncation radius for both lattices
 * * `layer1_vectors`   – flat [lv1x, lv1y, lv2x, lv2y] for layer 1
 * * `layer2_vectors`   – flat [lv1x, lv1y, lv2x, lv2y] for layer 2 (currently unused –
 *                        both lattices are generated from layer1_vectors, matching the JS)
 * * `tol`              – decimal precision for distance grouping
 * * `angle_start`      – lower bound of the angle search window (degrees, exclusive)
 * * `angle_end`        – upper bound of the angle search window (degrees, exclusive)
 *
 * # Returns
 * A flat `Vec<f64>` with **7 values per result row**, in this order:
 *   `angle_deg | angle_rad | ll1 | ll2 | ul1 | ul2 | cells`
 * Total length = 7 × N results.
 */
export function find_values(radius: number, layer1_vectors: Float64Array, _layer2_vectors: Float64Array, tol: number, angle_start: number, angle_end: number): Float64Array;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly find_values: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
