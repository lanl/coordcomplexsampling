import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import argparse
import pathlib


def perform_swap(swapdict, writeout=False, output_path="output_sdfs",
                 call=0, skip_checks=False):
    """perform the swap

    Args:
        swapdict (dict): swap dictionary
        writeout (bool, optional): write sdf file. Defaults to False.
        output_path (str, optional): where to write out. Defaults to "output_sdfs".
        call (int, optional): number of call attempt. Defaults to 0.
        skip_checks (bool, optional): return regardless of distances.

    Returns:
        sdfstr : init/final functionalized sdf str.
    """
    import os

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"

    from architector import convert_io_molecule
    import pathlib

    init_mol = convert_io_molecule(swapdict["init_sdf"])
    init_mol.charge = 0
    init_mol.uhf = 0
    final_mol = convert_io_molecule(swapdict["final_sdf"])
    final_mol.charge = 0
    final_mol.uhf = 0
    # Swap metal(s)
    for i, metind in enumerate(swapdict["metal_inds"]):
        init_mol.ase_atoms[metind].symbol = swapdict["metal_syms"][i]
        init_mol.atom_types[metind] = swapdict["metal_syms"][i]
        final_mol.ase_atoms[metind].symbol = swapdict["metal_syms"][i]
        final_mol.atom_types[metind] = swapdict["metal_syms"][i]
    # Functionalize 3D
    if swapdict["swap_remove_inds"] is not None:
        init_mol.functionalize_3D(
            functional_groups=swapdict["swap_functional_groups"],
            functionalization_inds=swapdict["swap_functionalization_inds"],
            functional_group_mol_inds=swapdict[
                "swap_functional_group_mol_inds"
            ],
            remove_inds=swapdict["swap_remove_inds"],
            bond_orders=swapdict["swap_bond_orders"],
            remove_hydrogens_when_adding=swapdict[
                "swap_remove_hydrogens_when_adding"
            ],
            xtb_opt=swapdict["swap_xtb_opt"],
        )
        final_mol.functionalize_3D(
            functional_groups=swapdict["swap_functional_groups"],
            functionalization_inds=swapdict["swap_functionalization_inds"],
            functional_group_mol_inds=swapdict[
                "swap_functional_group_mol_inds"
            ],
            remove_inds=swapdict["swap_remove_inds"],
            bond_orders=swapdict["swap_bond_orders"],
            remove_hydrogens_when_adding=swapdict[
                "swap_remove_hydrogens_when_adding"
            ],
            xtb_opt=swapdict["swap_xtb_opt"],
        )
    init_mol.dist_sanity_checks(
        smallest_dist_cutoff=0.7,
        min_dist_cutoff=10,
    )  # Check that atoms are not overlapping
    init_mol.graph_sanity_checks()
    if init_mol.dists_sane or skip_checks:
        init_sdf = init_mol.write_sdf("init", writestring=True)
        final_sdf = final_mol.write_sdf("final", writestring=True)
        outstr = init_sdf + final_sdf
        outname = swapdict["rxn_uid"] + ".sdf"
        if writeout:
            tp = pathlib.Path(output_path)
            if not tp.exists():
                tp.mkdir(exist_ok=True, parents=True)
            with open(tp / outname, "w") as file1:
                file1.write(outstr)
            return "Done"
        else:
            # print(outname)
            return outstr
    elif call < 11:
        out = perform_swap(
            swapdict, writeout=writeout, output_path=output_path, call=call + 1
        )
        return out
    else:
        if writeout:
            outname = swapdict["rxn_uid"] + ".txt"
            tp = pathlib.Path(output_path)
            if not tp.exists():
                tp.mkdir(exist_ok=True, parents=True)
            with open(tp / outname, "w") as file1:
                file1.write(
                    "RXN: "
                    + swapdict["rxn_uid"]
                    + " FAILED AFTER 10X attempts."
                )
            return "Done"
        else:
            return None


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample_dataframe", default="to_gen_structures_parallel.pkl"
    )
    parser.add_argument("--output_path", default="output_sdfs")
    parser.add_argument("--nprocs", default=12, type=int)
    parser.add_argument("--chunksize", default=5000, type=int)
    parser.add_argument("--writeout", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    outpath = pathlib.Path(args.output_path)

    print(args.__dict__)

    if args.writeout:
        done_lst = [x.name.replace(".sdf", "") for x in outpath.rglob("*.sdf")]
        done_lst += [x.name.replace(".txt", "") for x in outpath.rglob("*.txt")]
    else:
        done_lst = []

    df = pd.read_pickle(args.sample_dataframe)

    with tqdm(total=(df.shape[0] - len(done_lst))) as pbar:
        with ProcessPoolExecutor(max_workers=args.nprocs) as exe:
            futs = []
            donedicts = []
            for i, row in df.iterrows():
                if row["rxn_uid"] not in done_lst:
                    indict = row.to_dict()
                    # swapdict, writeout=False, output_path='output_sdfs'
                    fut = exe.submit(
                        perform_swap,
                        indict,
                        writeout=args.writeout,
                        output_path=args.output_path,
                    )
                    fut.savename = row["rxn_uid"] + ".sdf"
                    futs.append(fut)
            for x in as_completed(futs):
                if not args.writeout:
                    donedicts.append(
                        {"sdf_name": x.savename, "sdf_file": x.result()}
                    )
                pbar.update(1)

    if len(donedicts) > 0:
        outdf = pd.DataFrame(donedicts)
        outdf.to_pickle("all_combined.pkl")
