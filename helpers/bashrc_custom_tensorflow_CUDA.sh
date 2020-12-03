#!/bin/sh

    ## aliroot env
    export SINGULARITY_PREFIX=/home/marsland/containers/tensorflow_latest_gpu_CUDA/
    export SINGULARITY_CACHEDIR=/home/marsland/containers/tensorflow_latest_gpu_CUDA/CACHEDIR
    cd ${SINGULARITY_PREFIX}/tensorflow_latest-gpuCUDA/alicesw/
    export ALIBUILD_WORK_DIR="${SINGULARITY_PREFIX}/tensorflow_latest-gpuCUDA/alicesw/sw"
    eval "`alienv shell-helper`"
    #alias aliload6='alienv load AliPhysics/latest-master-next-root6; alienv list'
    alias aliload6='alienv load AliPhysics/latest-aliroot6-user-next-root6; alienv list'
    alias aliload5='alienv load AliPhysics/latest-aliroot5-user; alienv list'

    ###########################################################################################################
    ############################################# FROM bashrc #################################################
    ###########################################################################################################

    export TIdentityDIR=/home/marsland/Desktop/ubuntu_desktop/workdir/TIdentity_ND
    export RUN_ON_GRID_DIR=/home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID
    export NOTES=/home/marsland/Desktop/ubuntu_desktop/github/alice-tpc-notes
    export NOTESData=/home/marsland/Desktop/ubuntu_desktop/workdir/NOTESData


    ## some more ls aliases
    alias ll='ls -l'
    alias lh='ls -ltrh'
    alias lhs='ls -ltrhS'
    alias la='ls -A'
    alias l='ls -CF'
    alias lsdots='ls -ald .*'
    alias less='less -R'

    ## gsi related
    alias sshgsi='ssh -XY marsland@lxpool.gsi.de'
    alias sshcern='ssh -XY marsland@lxplus.cern.ch'
    alias sshhei='ssh -XY arslandok@pi0.physi.uni-heidelberg.de'
    alias umountgsi=". /home/marsland/bin/umountGSI.sh"
    alias mountgsi=". /home/marsland/bin/mountDisks.sh"
    alias cdafs="cd /afs/cern.ch/user/m/marsland"
    alias cdn="cd /lustre/nyx/alice/users/marsland"
    alias cdwork='cd /home/marsland/Desktop/ubuntu_desktop/workdir'

    ## don't put duplicate lines or lines starting with space in the history.
    HISTCONTROL=ignoreboth

    ## append to the history file, don't overwrite it
    shopt -s histappend

    ## for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
    HISTSIZE=1000
    HISTFILESIZE=2000

    ## check the window size after each command and, if necessary, update the values of LINES and COLUMNS.
    shopt -s checkwinsize

    ## make less more friendly for non-text input files, see lesspipe(1)
    [ -x /usr/bin/lesspipe ] && eval "$(SHELL=/bin/sh lesspipe)"

    ## set variable identifying the chroot you work in (used in the prompt below)
    if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
        debian_chroot=$(cat /etc/debian_chroot)
    fi

    ## set a fancy prompt (non-color, unless we know we "want" color)
    case "$TERM" in
        xterm-color|*-256color) color_prompt=yes;;
    esac

    if [ -n "$force_color_prompt" ]; then
        if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
      # We have color support; assume it's compliant with Ecma-48
      # (ISO/IEC-6429). (Lack of such support is extremely rare, and such
      # a case would tend to support setf rather than setaf.)
      color_prompt=yes
        else
      color_prompt=
        fi
    fi

    if [ "$color_prompt" = yes ]; then
        PS1='<aliSing> ${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
    else
        PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
    fi

    unset color_prompt force_color_prompt

    # If this is an xterm set the title to user@host:dir
    case "$TERM" in
    xterm*|rxvt*)
        PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
        ;;
    *)
        ;;
    esac

    # enable color support of ls and also add handy aliases
    if [ -x /usr/bin/dircolors ]; then
        test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
        alias ls='ls --color=auto'
        alias grep='grep --color=auto'
        alias fgrep='fgrep --color=auto'
        alias egrep='egrep --color=auto'
    fi

    # colored GCC warnings and errors
    export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'
