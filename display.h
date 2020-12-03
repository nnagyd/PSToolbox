#ifndef DISPLAY_HEADER
#define DISPLAY_HEADER

/*LWP*/

//void LWP::Plot_p() {
//	printf("\n Plotting...");
//	int ll = data.size();
//	vector<double> t(ll), pe(ll), pv(ll);
//	for (int i = 0; i < ll; i++) {
//		t.at(i) = data.at(i).at(0);
//		pe.at(i) = data.at(i).at(1) / 1.e5;
//		pv.at(i) = data.at(i).at(2) / 1.e5;
//		/*printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));*/
//	}
//	plt::plot(t, pe, "k", t, pv, "b");
//	plt::xlabel("t (s)");
//	plt::ylabel("p (bar)");
//	plt::title(name);
//	plt::grid(true);
//	plt::show();
//}
//
//void LWP::Plot_v() {
//	printf("\n Plotting...");
//	int ll = data.size();
//	vector<double> t(ll), xe(ll), xv(ll);
//	for (int i = 0; i < ll; i++) {
//		t.at(i) = data.at(i).at(0);
//		xe.at(i) = data.at(i).at(3);
//		xv.at(i) = data.at(i).at(4);
//		/*printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));*/
//	}
//	plt::plot(t, xe, "k", t, xv, "b");
//	plt::xlabel("t (s)");
//	plt::ylabel("v (m/s)");
//	plt::title(name);
//	plt::grid(true);
//	plt::show();
//}
//
//void LWP::Plot_T() {
//	printf("\n Plotting...");
//	int ll = data.size();
//	vector<double> t(ll), xe(ll), xv(ll);
//	for (int i = 0; i < ll; i++) {
//		t.at(i) = data.at(i).at(0);
//		xe.at(i) = data.at(i).at(5) -273.15;
//		xv.at(i) = data.at(i).at(6) -273.15;
//		/*printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));*/
//	}
//	plt::plot(t, xe, "k", t, xv, "b");
//	plt::xlabel("t (s)");
//	plt::ylabel("T (C)");
//	plt::title(name);
//	plt::grid(true);
//	plt::show();
//}
//
//void LWP::Plot_rho() {
//	printf("\n Plotting...");
//	int ll = data.size();
//	vector<double> t(ll), xe(ll), xv(ll);
//	for (int i = 0; i < ll; i++) {
//		t.at(i) = data.at(i).at(0);
//		xe.at(i) = data.at(i).at(7) / 1.e5;
//		xv.at(i) = data.at(i).at(8) / 1.e5;
//		/*printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));*/
//	}
//	plt::plot(t, xe, "k", t, xv, "b");
//	plt::xlabel("t (s)");
//	plt::ylabel("rho (kg/m3)");
//	plt::title(name);
//	plt::grid(true);
//	plt::show();
//}
//
//void LWP::Plot_mp() {
//	printf("\n Plotting...");
//	int ll = data.size();
//	vector<double> t(ll), xe(ll), xv(ll);
//	for (int i = 0; i < ll; i++) {
//		t.at(i) = data.at(i).at(0);
//		xe.at(i) = data.at(i).at(9) / 1.e5;
//		xv.at(i) = data.at(i).at(10) / 1.e5;
//		/*printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));*/
//	}
//	plt::plot(t, xe, "k", t, xv, "b");
//	plt::xlabel("t (s)");
//	plt::ylabel("mp (kg/s)");
//	plt::title(name);
//	plt::grid(true);
//	plt::show();
//}



#endif // !DISPLAY_HEADER
