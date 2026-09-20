#!/usr/bin/env bash
# Source this file:
#   source slurm_workdir_tools.sh
#
# These helpers operate by WorkDir, which is safer than matching job names.

q_by_root () {
    local root="$1"
    squeue -u "$USER" -h -O JobID:20,Name:25,State:15,Reason:25,WorkDir:250 \
    | awk -v root="$root" '$5 ~ "^"root {print}'
}

cancel_by_root () {
    local root="$1"
    squeue -u "$USER" -h -O JobID:20,WorkDir:250 \
    | awk -v root="$root" '$2 ~ "^"root {print $1}' \
    | xargs -r scancel
}

cancel_by_dir () {
    local dir
    dir="$(realpath "$1")"
    squeue -u "$USER" -h -O JobID:20,WorkDir:250 \
    | awk -v dir="$dir" '$2 == dir {print $1}' \
    | xargs -r scancel
}

dupes_by_root () {
    local root="$1"
    squeue -u "$USER" -h -O JobID:20,Name:25,State:15,Reason:25,WorkDir:250 \
    | awk -v root="$root" '$5 ~ "^"root {print $5}' \
    | sort | uniq -c | awk '$1 > 1'
}

submit_clean_dirs () {
    # Usage:
    #   submit_clean_dirs "/path/to/root" "level*_*-ads" "first_vasp.run"
    local root="$1"
    local pattern="${2:-level*_*-ads}"
    local submit="${3:-first_vasp.run}"

    cd "$root" || return 1

    for d in $pattern; do
        [ -d "$d" ] || continue
        [ -f "$d/DO_NOT_SUBMIT_OVERLAP" ] && {
            echo "SKIP quarantined: $d"
            continue
        }
        [ -f "$d/$submit" ] || {
            echo "SKIP missing submit file: $d/$submit"
            continue
        }

        echo "Submitting $d"
        (
            cd "$d" || exit 1
            rm -f WAVECAR CHGCAR
            sbatch "$submit"
        )
    done
}
