% error computation

IDconv1Simulink=interp1(out.tout,out.xout{15}.Values.Data(:,1),conv1Simout.tsim);
IQconv1Simulink=interp1(out.tout,out.xout{15}.Values.Data(:,2),conv1Simout.tsim);

VDconv1Simulink=interp1(out.tout,networkData(5,:),conv1Simout.tsim);
VQconv1Simulink=interp1(out.tout,networkData(6,:),conv1Simout.tsim);

error.normalizedMaxDeviationIDconv1=max(abs(conv1Simout.x(:,1)-IDconv1Simulink))/(max(IDconv1Simulink)-min(IDconv1Simulink));
error.normalizedRMSdeviationIDconv1=sqrt(sum((IDconv1Simulink-[conv1Simout.x(:,1)]).^2)/length(conv1Simout.tsim))/(max(IDconv1Simulink)-min(IDconv1Simulink));
error.meanAbsoluteErrorIDconv1=sum(abs(conv1Simout.x(:,1)-IDconv1Simulink))/length(conv1Simout.x(:,1));

error.normalizedMaxDeviationIQconv1=max(abs(conv1Simout.x(:,2)-IQconv1Simulink))/(max(IQconv1Simulink)-min(IQconv1Simulink));
error.normalizedRMSdeviationIQconv1=sqrt(sum((IDconv1Simulink-[conv1Simout.x(:,2)]).^2)/length(conv1Simout.tsim))/(max(IQconv1Simulink)-min(IQconv1Simulink));
error.meanAbsoluteErrorIQconv1=sum(abs(conv1Simout.x(:,2)-IQconv1Simulink))/length(conv1Simout.x(:,2));

error.normalizedMaxDeviationVDconv1=max(abs(conv1Simout.y(:,1)-VDconv1Simulink))/(max(VDconv1Simulink)-min(VDconv1Simulink));
error.normalizedRMSdeviationVDconv1=sqrt(sum((VDconv1Simulink-[conv1Simout.y(:,1)]).^2)/length(conv1Simout.tsim))/(max(VDconv1Simulink)-min(VDconv1Simulink));
error.meanAbsoluteErrorVDconv1=sum(abs(conv1Simout.y(:,1)-VDconv1Simulink))/length(conv1Simout.y(:,1));

error.normalizedMaxDeviationVQconv1=max(abs(conv1Simout.y(:,2)-VQconv1Simulink))/(max(VQconv1Simulink)-min(VQconv1Simulink));
error.normalizedRMSdeviationVQconv1=sqrt(sum((VQconv1Simulink-[conv1Simout.y(:,2)]).^2)/length(conv1Simout.tsim))/(max(VQconv1Simulink)-min(VQconv1Simulink));
error.meanAbsoluteErrorVQconv1=sum(abs(conv1Simout.y(:,2)-VQconv1Simulink))/length(conv1Simout.y(:,2));


 errorMetrics(k,:)=[timeLinearization,error.meanAbsoluteErrorIDconv1,error.normalizedMaxDeviationIDconv1,error.meanAbsoluteErrorIQconv1,error.normalizedMaxDeviationIQconv1,error.meanAbsoluteErrorVDconv1,error.normalizedMaxDeviationVDconv1,error.meanAbsoluteErrorVQconv1,error.normalizedMaxDeviationVQconv1];

 errorTable=table(errorMetrics(:,1),errorMetrics(:,2),errorMetrics(:,3),errorMetrics(:,4),errorMetrics(:,5),errorMetrics(:,6),errorMetrics(:,7),errorMetrics(:,8),errorMetrics(:,9),'VariableName',{'Time initialization','IDconv1 MAE','IDconv1 norm. max. E','IQconv1 MAE','IQconv1 norm. max. E','VDconv1 MAE','VDconv1 norm. max. E','VQconv1 MAE','VQconv1 norm. max. E'});
