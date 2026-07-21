#!/bin/bash -ue
echo "step 1, sample [[aT5KCaxc:[525025134430777892], vohFRHsV:[288073750975865022], vNAT43RD:[939339087586825539], mIlFiO6w:[998115530373274444], Ip5yo1fu:[706496068017274183], xdeUzdZ4:[727968464493146190], CS92or1c:[720626292830168825], FILES:[[/msm_home/runs/XmwMJQzm/_metasmith/task/data/Ymf7qCbcMqFp/scadc_main.txt], [/home/tony/agentic_workspace/projects/fabfos/scadc/.awm/data/ref/derived/mnxref-4_5/solve_directed/reff_axes_report.tsv], [/home/tony/agentic_workspace/projects/fabfos/scadc/.awm/data/ref/derived/mnxref-4_5/solve_directed/ieff_axes_report.tsv], [/home/tony/agentic_workspace/data/scadc/fabfos_2026/orfs_199.faa], [/home/tony/agentic_workspace/projects/fabfos/scadc/methods/.runs/scadc_main/refs/ecspr_frozen_null], [/msm_home/runs/XmwMJQzm/_metasmith/task/data/IvyoACq9xkYX/ecspr.oci], [/msm_home/runs/XmwMJQzm/_metasmith/task/data/GWVC6WPW1X1W/ecspr_significance.py]]]]"
echo "significance"
echo "res 4/6 GB/1" >>.command.metadata
echo "lin [{\"aT5KCaxc\":[525025134430777892],\"vohFRHsV\":[288073750975865022],\"vNAT43RD\":[939339087586825539],\"mIlFiO6w\":[998115530373274444],\"Ip5yo1fu\":[706496068017274183],\"xdeUzdZ4\":[727968464493146190],\"CS92or1c\":[720626292830168825],\"FILES\":[[\"/msm_home/runs/XmwMJQzm/_metasmith/task/data/Ymf7qCbcMqFp/scadc_main.txt\"],[\"/home/tony/agentic_workspace/projects/fabfos/scadc/.awm/data/ref/derived/mnxref-4_5/solve_directed/reff_axes_report.tsv\"],[\"/home/tony/agentic_workspace/projects/fabfos/scadc/.awm/data/ref/derived/mnxref-4_5/solve_directed/ieff_axes_report.tsv\"],[\"/home/tony/agentic_workspace/data/scadc/fabfos_2026/orfs_199.faa\"],[\"/home/tony/agentic_workspace/projects/fabfos/scadc/methods/.runs/scadc_main/refs/ecspr_frozen_null\"],[\"/msm_home/runs/XmwMJQzm/_metasmith/task/data/IvyoACq9xkYX/ecspr.oci\"],[\"/msm_home/runs/XmwMJQzm/_metasmith/task/data/GWVC6WPW1X1W/ecspr_significance.py\"]]}]" >>.command.metadata
echo "fmt 2" >>.command.metadata
cat /home/tony/agentic_workspace/projects/fabfos/scadc/methods/.runs/scadc_main/agent_home/runs/XmwMJQzm/workflow.step_1.meta >>.command.metadata
echo "inp aT5KCaxc,vohFRHsV,vNAT43RD,mIlFiO6w,Ip5yo1fu,xdeUzdZ4,CS92or1c" >>.command.metadata
echo "out pKixB8Yw,i3oS98wX" >>.command.metadata
b1="/home/tony/agentic_workspace"
echo "--bind $b1:$b1" >.command.binds

CONTAINER=/msm_home
DIRECT=/home/tony/agentic_workspace/projects/fabfos/scadc/methods/.runs/scadc_main/agent_home
function bootstrap {
	if [ -e $CONTAINER ]; then
		$CONTAINER/lib/msm_bootstrap $@
	elif [ -e $DIRECT ]; then
		$DIRECT/lib/msm_bootstrap $@
	else
		echo "critical error: could not find metasmith bootstrap script"
	fi
}

bootstrap /home/tony/agentic_workspace/projects/fabfos/scadc/methods/.runs/scadc_main/agent_home/runs/XmwMJQzm "1" Cosmos
[ -e .command.success ] && exit 0 || exit 1
