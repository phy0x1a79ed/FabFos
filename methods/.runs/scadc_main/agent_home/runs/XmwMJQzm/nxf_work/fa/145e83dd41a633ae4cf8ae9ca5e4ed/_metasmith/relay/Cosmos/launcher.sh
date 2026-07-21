SCRIPT=$1
PIDF=$2
DONEF=$3
cd /home/tony/agentic_workspace/projects/fabfos/scadc/methods/.runs/scadc_main/agent_home/runs/XmwMJQzm/nxf_work/fa/145e83dd41a633ae4cf8ae9ca5e4ed
bash /home/tony/agentic_workspace/projects/fabfos/scadc/methods/.runs/scadc_main/agent_home/runs/XmwMJQzm/nxf_work/fa/145e83dd41a633ae4cf8ae9ca5e4ed/_metasmith/relay/Cosmos/$SCRIPT &
cd /home/tony/agentic_workspace/projects/fabfos/scadc/methods/.runs/scadc_main/agent_home/runs/XmwMJQzm/nxf_work/fa/145e83dd41a633ae4cf8ae9ca5e4ed/_metasmith/relay/Cosmos
PID=$!
echo $PID > $PIDF
wait $PID
STATUS=$?
rm $PIDF
rm $SCRIPT
echo $STATUS > $DONEF