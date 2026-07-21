params.testSpread=1
params.home = '/home/tony/agentic_workspace/projects/fabfos/scadc/methods/.runs/scadc_main/agent_home'
params.workspace = "${params.home}/runs/XmwMJQzm"
params.bootstrap_def = '''
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
'''

def in(f, l) {
    def rows = Channel.fromPath(f).splitCsv(header: false)
    if (f in l) {
        rows = Channel.fromList(l[f]).merge(rows)
    }
    return rows.map { row ->
        if (row.size()>1) {
            def (ri, rx) = row
            return tuple(ri, file(rx))
        } else {
            def i = [:]
            return tuple(i, file(row[0]))
        }
    }
}


process p01__significance {
	label 'xD4CvuFyIx'
input:
	tuple val(index),path(_01),path(_02),path(_03),path(_04),path(_05),path(_06),path(_07)
output:
	tuple val(index),path("*-1.*-pKixB8Yw.tsv")
	tuple val(index),path("*-1.*-i3oS98wX.tsv")
script:
"""
echo "step 1, sample $index"
echo "significance"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "fmt 2" >>.command.metadata
cat ${params.workspace}/workflow.step_1.meta >>.command.metadata
echo "inp aT5KCaxc,vohFRHsV,vNAT43RD,mIlFiO6w,Ip5yo1fu,xdeUzdZ4,CS92or1c" >>.command.metadata
echo "out pKixB8Yw,i3oS98wX" >>.command.metadata
b1="/home/tony/agentic_workspace"
echo "--bind \$b1:\$b1" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "1" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries { k, v -> [k, v.sort()] }}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-pKixB8Yw.tsv" "1-1-1.test$hash-i3oS98wX.tsv"
"""
}

workflow {
main:
o = new Orchestrator(Channel.fromList([null])) // cant create channels in groovy
_lf = new groovy.json.JsonSlurper().parseText(file("workflow.lineage_of_given.json").text)
l = _lf.lineage
o.seedParents(_lf.child2parent)
_xdeUzdZ4 = (o.postIn([in("inputs/xdeUzdZ4", l)], ["xdeUzdZ4"]))[0] // containers::ecspr.oci
_Ip5yo1fu = (o.postIn([in("inputs/Ip5yo1fu", l)], ["Ip5yo1fu"]))[0] // ecspr::frozen_null
_aT5KCaxc = (o.postIn([in("inputs/aT5KCaxc", l)], ["aT5KCaxc"]))[0] // fosmids::recovery_experiment
_CS92or1c = (o.postIn([in("inputs/CS92or1c", l)], ["CS92or1c"]))[0] // lib::ecspr_significance.py
_vNAT43RD = (o.postIn([in("inputs/vNAT43RD", l)], ["vNAT43RD"]))[0] // ecspr::ieff_axes_report
_vohFRHsV = (o.postIn([in("inputs/vohFRHsV", l)], ["vohFRHsV"]))[0] // ecspr::reff_axes_report
_mIlFiO6w = (o.postIn([in("inputs/mIlFiO6w", l)], ["mIlFiO6w"]))[0] // sequences::open_reading_frames
k = ['pKixB8Yw', 'i3oS98wX']
(_pKixB8Yw, _i3oS98wX) = o.post(o.asStreams(p01__significance(o.group('aT5KCaxc', [_aT5KCaxc, _vohFRHsV, _vNAT43RD, _mIlFiO6w, _Ip5yo1fu, _xdeUzdZ4, _CS92or1c], k, 1))), k)

publish:
_i3oS98wX = o.publish(_i3oS98wX)
_pKixB8Yw = o.publish(_pKixB8Yw)
}

output {
	_pKixB8Yw{
		path 'ecspr-reff_significance'
		index { path '_manifests/ecspr-reff_significance.pKixB8Yw.K5qZE5befF.json' }
	}
	_i3oS98wX{
		path 'ecspr-ieff_significance'
		index { path '_manifests/ecspr-ieff_significance.i3oS98wX.wKt2UY8OOK.json' }
	}
}