SIM_EXP_NAME = test-sim-6-8-1
awsbucket = s3://cytof-repfam/$(SIM_EXP_NAME)/viz

sendimg:
	aws s3 sync viz $(awsbucket)

getimg:
	aws s3 sync $(awsbucket) viz
