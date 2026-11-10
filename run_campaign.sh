#!/usr/bin/env bash
# Copy a parameter file, apply user inputs, and run each case in its own folder.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$ROOT/parameters"
SAVE="$ROOT/campaign_runs/$(date +%Y%m%d-%H%M%S)"
MODE=interactive FILE= NAME=current MODEL= JOBS=1 QUIET=0 DRY=0
CASES=() OVERRIDES=() DIRS=() MODES=()

fail() { echo "ERROR: $*" >&2; exit 2; }
trim() {
    local s="$1"
    s="${s#"${s%%[![:space:]]*}"}"
    s="${s%"${s##*[![:space:]]}"}"
    printf '%s' "$s"
}
usage() {
    cat <<'HELP'
STRATA-NOVA campaign runner

  ./run_campaign.sh                         Enter parameter=value changes
  ./run_campaign.sh --current               Run the current parameters
  ./run_campaign.sh --set csm=none --set simulation_mode=single
  ./run_campaign.sh --file example/campaign.example --jobs 4

Options:
  --model NAME    Select both bundled DDT model files, e.g. n20 or n1600c
  --list-models   List the available bundled DDT models
  --base FILE     Start from another complete parameter file
  --set KEY=VALUE Override any parameter in the base or root template (repeatable)
  --name NAME     Name a single case (default: current)
  --file FILE     Cases: name | key=value | key=value (one case per line)
  --current       Run once without interactive prompts
  --save DIR      Results folder (default: campaign_runs/<timestamp>)
  --jobs N        Concurrent cases (default: 1; solver logs stay in run.log)
  --quiet         Keep solver output in run.log
  --dry-run       Prepare and check inputs without building or running
  -h, --help      Show this help

String values may be quoted or unquoted. Arrays use commas, e.g. 2, 7, 12.
Model paths are relative to the repository root, or absolute.
--set changes apply to every case; case-file changes take precedence.
Each case keeps parameters.used, run.log, exit_status and output/.
The runner chooses output paths and never edits your base parameter file.
A dry run checks input syntax/keys; the solver validates physical values.
HELP
}

# Read an active assignment. Strip comments only outside quoted strings.
get_value() {
    awk -v key="$2" '
      /^[[:space:]]*[#!]/ {next}
      {
        eq=index($0,"="); name=substr($0,1,eq-1); gsub(/[ \t]/,"",name)
        if (eq && name==key) {
          value=substr($0,eq+1); quoted=0
          for(i=1;i<=length(value);i++) {
            c=substr(value,i,1)
            if(c=="\"") quoted=!quoted
            if(!quoted && (c=="!" || c=="#")) {value=substr(value,1,i-1); break}
          }
          gsub(/^[ \t]+|[ \t]+$/,"",value)
          gsub(/^"|"$/,"",value); print value; exit
        }
      }' "$1"
}

# Apply every override in one pass, using the templates as the key/type list.
prepare_parameters() {
    local file="$1" change temporary
    shift
    for change in "$@"; do
        [[ "$change" == *[![:space:]]* ]] || fail "Empty parameter override"
        [[ "$change" != *$'\n'* ]] || fail "Each override must occupy one line"
    done
    temporary="$(mktemp "$file.XXXXXX")"
    if ! CHANGES="$(printf '%s\n' "$@")" awk '
      function trim(s) {gsub(/^[ \t]+|[ \t]+$/,"",s); return s}
      function error(message) {print "ERROR: " message > "/dev/stderr"; exit 2}
      {
        line=$0; sub(/^[ \t]*[#!]?[ \t]*/,"",line)
        eq=index(line,"="); key=trim(substr(line,1,eq-1))
        if(eq && key ~ /^[A-Za-z_][A-Za-z0-9_]*$/) {
          known[key]=1
          if(substr(line,eq+1) ~ /^[ \t]*"/) string_key[key]=1
        }
        if(NR==FNR) next
        lines[++count]=$0
      }
      END {
        n=split(ENVIRON["CHANGES"],changes,"\n")
        for(i=1;i<=n;i++) {
          line=trim(changes[i]); eq=index(line,"=")
          key=trim(substr(line,1,eq-1)); value=trim(substr(line,eq+1))
          if(!eq || key !~ /^[A-Za-z_][A-Za-z0-9_]*$/ || value=="") error("Expected parameter=value: " line)
          if(key=="outdir" || key=="staged_workdir") error("Use --save to choose output folders.")
          if(!(key in known)) error("Unknown parameter " key ". Use ? to see the template.")
          if(string_key[key] && substr(value,1,1)!="\"") value="\"" value "\""
          if(!(key in values)) order[++changed]=key
          values[key]=value
        }
        for(i=1;i<=count;i++) {
          eq=index(lines[i],"="); key=trim(substr(lines[i],1,eq-1))
          if(lines[i] !~ /^[ \t]*[#!]/ && eq && key in values) {
            if(!written[key]++) print key " = " values[key]
          } else print lines[i]
        }
        for(i=1;i<=changed;i++) {
          key=order[i]
          if(!written[key]) print key " = " values[key]
        }
      }' "$ROOT/parameters" "$BASE" > "$temporary"; then
        rm -f "$temporary"
        return 2
    fi
    mv "$temporary" "$file"
}

list_models() {
    local profile model
    for profile in "$ROOT"/profiles/ddt_cases/wdprofile_ddt_2013_*.dat; do
        model="${profile##*ddt_2013_}"; model="${model%.dat}"
        [[ ! -f "$profile" || ! -f "$ROOT/profiles/ddt_cases/wdcomp_ddt_2013_$model.dat" ]] || printf '%s\n' "$model"
    done
}

while (($#)); do
    case "$1" in
      --base|--save|--name|--file|--jobs|--set|--model)
        (($# >= 2)) || fail "$1 needs a value"
        case "$1" in
          --base) BASE="$2" ;; --save) SAVE="$2" ;; --name) NAME="$2" ;;
          --file) [[ "$MODE" != current ]] || fail "Choose --file or --current"; MODE=file; FILE="$2" ;;
          --jobs) JOBS="$2" ;; --model) MODEL="${2,,}" ;;
          --set) OVERRIDES+=("$2") ;;
        esac
        shift 2 ;;
      --list-models) list_models; exit 0 ;;
      --current) [[ "$MODE" != file ]] || fail "Choose --file or --current"; MODE=current; shift ;;
      --quiet) QUIET=1; shift ;; --dry-run) DRY=1; shift ;;
      -h|--help) usage; exit 0 ;;
      *) fail "Unknown option '$1'; see --help" ;;
    esac
 done
if [[ -n "$MODEL" ]]; then
    [[ "$MODEL" =~ ^n[0-9]+[a-z]*$ && -f "$ROOT/profiles/ddt_cases/wdprofile_ddt_2013_$MODEL.dat" && -f "$ROOT/profiles/ddt_cases/wdcomp_ddt_2013_$MODEL.dat" ]] || fail "Unknown bundled model '$MODEL'. Use --list-models."
    OVERRIDES=("profile_name=profiles/ddt_cases/wdprofile_ddt_2013_$MODEL.dat"
        "comp_profile_name=profiles/ddt_cases/wdcomp_ddt_2013_$MODEL.dat" "${OVERRIDES[@]}")
fi
if [[ "$MODE" == interactive && ${#OVERRIDES[@]} -gt 0 ]]; then MODE=current; fi
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || fail "--jobs must be a positive integer"
[[ -f "$BASE" ]] || fail "Base parameter file not found: $BASE"
BASE="$(cd "$(dirname "$BASE")" && pwd)/$(basename "$BASE")"

if [[ "$MODE" == interactive ]]; then
    [[ -t 0 ]] || fail "Use --current, --set, or --file for non-interactive runs"
    echo "Starting from $BASE"
    echo "Enter parameter=value, ? to view the template, or a blank line to finish."
    echo "Examples (one per line): csm=active, mass_scheme=constant_density, csm_mass_msun=0.01"
    echo "For a new model, set both profile_name and comp_profile_name."
    while read -r -p '> ' line; do
        line="$(trim "$line")"
        [[ -n "$line" ]] || break
        if [[ "$line" == '?' ]]; then cat "$BASE"; else OVERRIDES+=("$line"); fi
    done
    read -r -p "Case name [$NAME]: " line || fail "Input closed"
    [[ -z "$line" ]] || NAME="$line"
    read -r -p "Results folder [$SAVE]: " line || fail "Input closed"
    [[ -z "$line" ]] || SAVE="$line"
fi
if [[ "$MODE" == file ]]; then
    [[ -f "$FILE" ]] || fail "Case file not found: $FILE"
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="$(trim "${line%$'\r'}")"
        [[ -z "$line" || "$line" == \#* ]] || CASES+=("$line")
    done < "$FILE"
else
    CASES=("$NAME")
fi
((${#CASES[@]})) || fail "No cases selected"
# Summaries include only cases selected in this invocation.
mkdir -p "$SAVE"
SAVE="$(cd "$SAVE" && pwd)"
[[ "$SAVE" != *\"* && "$SAVE" != *$'\n'* && "$SAVE" != *'|'* ]] || fail "Results path cannot contain quotes, newlines, or |"

# Prepare ALL inputs before any build or solver run; typos cannot start half a grid.
declare -A SEEN=()
for specification in "${CASES[@]}"; do
    IFS='|' read -r -a parts <<< "$specification"
    label="$(trim "${parts[0]}")"
    [[ "$label" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]*$ ]] || fail "Case names must start with a letter, digit or _ and use only letters, digits, _, . or -"
    dir="$SAVE/$label"
    [[ -z "${SEEN[$label]:-}" ]] || fail "Duplicate case name '$label' in this campaign"
    SEEN[$label]=1
    [[ ! -L "$dir" ]] || fail "Case folder cannot be a symbolic link: $dir"
    if [[ -d "$dir" ]]; then echo "Current folder is being overwritten: $dir"; fi
    mkdir -p "$dir"
    params="$dir/parameters.used"
    [[ "$BASE" != "$params" ]] || fail "This case would overwrite its base parameter file. Choose another --name or --save."
    prepare_parameters "$params" "${OVERRIDES[@]}" "${parts[@]:1}"
    mode="$(get_value "$params" simulation_mode)"
    mode="${mode:-single}"
    mode="${mode,,}"
    case "$mode" in
      single) key=outdir ;; staged_pipeline) key=staged_workdir ;;
      *) fail "Unknown simulation_mode '$mode' in '$label'" ;;
    esac
    # Output keys belong to the runner, and may be absent in a custom base.
    temporary="$(mktemp "$params.XXXXXX")"
    KEY="$key" VALUE="\"$dir/output\"" awk '
      BEGIN {key=ENVIRON["KEY"]}
      {name=$0; sub(/=.*/,"",name); gsub(/[ \t]/,"",name)
       if(name==key) next; print}
      END {print key " = " ENVIRON["VALUE"]}' "$params" > "$temporary"
    mv "$temporary" "$params"
    DIRS+=("$dir"); MODES+=("$mode")
    echo "Prepared $label ($mode): $dir"
done

if (( ! DRY )); then
    echo "Checking STRATA-NOVA build..."
    make --no-print-directory -s -C "$ROOT"
fi
(( JOBS == 1 )) || QUIET=1
run_case() {
    local dir="$1" mode="$2" status=0 result="$1/output"
    if (( DRY )); then
        echo 'dry-run: solver not started' > "$dir/run.log"
        echo dry-run > "$dir/exit_status"
        return 0
    fi
    mkdir -p "$result"
    local -a command=(env "OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}"
        "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}" "MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}"
        "SNEC_PARAMETER_FILE=$dir/parameters.used" "$ROOT/snec")
    if (( QUIET )); then
        (cd "$ROOT" && "${command[@]}") > "$dir/run.log" 2>&1 || status=$?
    else
        (cd "$ROOT" && "${command[@]}") 2>&1 | tee "$dir/run.log" || status=${PIPESTATUS[0]}
    fi
    [[ "$mode" != staged_pipeline ]] || result="$result/final"
    if (( status == 0 )) && [[ ! -s "$result/run_summary.dat" || ! -s "$result/lum_observed.dat" ]]; then
        echo "Missing completion products" >> "$dir/run.log"; status=1
    fi
    echo "$status" > "$dir/exit_status"
    echo "$(basename "$dir"): exit $status (log: $dir/run.log)"
    return "$status"
}

failures=0
for ((i=0; i<${#DIRS[@]}; )); do
    pids=()
    for ((n=0; n<JOBS && i<${#DIRS[@]}; n++, i++)); do
        if (( JOBS == 1 )); then
            run_case "${DIRS[i]}" "${MODES[i]}" || failures=$((failures+1))
        else
            run_case "${DIRS[i]}" "${MODES[i]}" & pids+=("$!")
        fi
    done
    for pid in "${pids[@]}"; do wait "$pid" || failures=$((failures+1)); done
done
{
    printf 'case\tdirectory\texit_status\tsimulation_mode\n'
    for ((i=0; i<${#DIRS[@]}; i++)); do
        label="$(basename "${DIRS[i]}")"
        printf '%s\t%s\t%s\t%s\n' "$label" "$label" "$(cat "${DIRS[i]}/exit_status")" "${MODES[i]}"
    done
} > "$SAVE/campaign_summary.tsv"
echo "Results: $SAVE (${#DIRS[@]} cases; $failures failed)"
(( failures == 0 ))
