#!/bin/bash
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# >>> agent setup commands
# <<<
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_DIR="./_metasmith/logs.$TIMESTAMP"
LOG_LATEST="./_metasmith/logs.latest"
mkdir -p $LOG_DIR
[ -e $LOG_LATEST ] && rm "$LOG_LATEST"; ln -s "./logs.$TIMESTAMP" "$LOG_LATEST"
[ -e workflow.params.yml ] || echo '{}' > workflow.params.yml
[ -e workflow.config.nf ] || touch workflow.config.nf
echo "start time was [$TIMESTAMP]"
export BINDS="--bind /home/tony/agentic_workspace:/home/tony/agentic_workspace"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
nohup ../../msm api run_workflow -a key=XmwMJQzm host=$(hostname) log_dir=$LOG_DIR stub_delay=${1:-0} >$LOG_DIR/agent.log 2>&1 &