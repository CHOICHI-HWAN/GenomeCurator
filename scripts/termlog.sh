# Keep a record of what was typed in a shell on this machine.
#
# Three machines run this cluster and work gets done by hand on all of them - a pull here, a worker restarted there - so
# when a machine turns out to be on the wrong code there is nothing to look back at. Sourcing this from ~/.bashrc or
# ~/.zshrc appends one line per command to a log, with the time, the machine, the directory and the exit status.
#
#   echo '. ~/BI_CLAUDE/GenomeCurator/scripts/termlog.sh' >> ~/.bashrc     # bash
#   echo '. ~/BI_CLAUDE/GenomeCurator/scripts/termlog.sh' >> ~/.zshrc      # zsh (macOS)
#
# The log is ~/.genomecurator/terminal-<host>.log, or wherever GC_TERMLOG points. It is this machine's own file: nothing
# is sent anywhere. Interactive shells only - a script that sources this is left alone.
#
# It records the command line as typed. A password or token typed on the command line would be recorded with it, so do
# not type secrets as arguments (and `git push` over HTTPS asks for them on a prompt, which is not recorded).

case $- in *i*) ;; *) return 0 ;; esac          # interactive shells only

: "${GC_TERMLOG:=$HOME/.genomecurator/terminal-$(hostname -s 2>/dev/null || hostname).log}"
mkdir -p "$(dirname "$GC_TERMLOG")" 2>/dev/null || return 0
[ -e "$GC_TERMLOG" ] || { : > "$GC_TERMLOG" && chmod 600 "$GC_TERMLOG"; } 2>/dev/null

__gc_write() {                                   # $1 exit status, $2 command line
    [ -n "$2" ] || return 0
    printf '%s\t%s\t%s\t%s\texit=%s\t%s\n' \
        "$(date '+%Y-%m-%dT%H:%M:%S')" "${HOSTNAME:-$(hostname -s 2>/dev/null)}" "${USER:-?}" "$PWD" "$1" "$2" \
        >> "$GC_TERMLOG" 2>/dev/null
}

if [ -n "$ZSH_VERSION" ]; then
    __gc_cmd=""
    __gc_preexec() { __gc_cmd="$1"; }
    __gc_precmd() {
        local st=$?
        [ -n "$__gc_cmd" ] && __gc_write "$st" "$__gc_cmd"
        __gc_cmd=""
    }
    autoload -Uz add-zsh-hook 2>/dev/null && { add-zsh-hook preexec __gc_preexec; add-zsh-hook precmd __gc_precmd; }
elif [ -n "$BASH_VERSION" ]; then
    # Read the last history entry at the prompt, after the command has run, so its exit status is the real one. The
    # history number tells a fresh command from a bare Enter, which would otherwise be logged again every prompt.
    # bash loads the history file after the startup files, so the number has to be taken at the first prompt rather than
    # here - otherwise that prompt logs the last command of an earlier session.
    __gc_n=""
    __gc_prompt() {
        local st=$? line n
        line=$(HISTTIMEFORMAT='' history 1) || return 0
        n=${line%%[!0-9 ]*}; n=${n// /}
        [ -n "$n" ] && [ "$n" != "$__gc_n" ] || return 0
        [ -z "$__gc_n" ] && { __gc_n="$n"; return 0; }      # the first prompt: nothing has been typed in this shell yet
        __gc_n="$n"
        __gc_write "$st" "${line#*[0-9]  }"
    }
    case ";$PROMPT_COMMAND;" in
        *";__gc_prompt;"*) ;;
        *) PROMPT_COMMAND="__gc_prompt${PROMPT_COMMAND:+;$PROMPT_COMMAND}" ;;
    esac
fi
