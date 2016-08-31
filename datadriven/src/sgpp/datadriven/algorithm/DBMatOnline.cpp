// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/base/exception/application_exception.hpp>

DBMatOnline::DBMatOnline (): offlineObject_(NULL){}

DBMatOnline::DBMatOnline(DBMatOffline* o) {
  readOffline(o);
}

DBMatOnline::~DBMatOnline(){
}

void DBMatOnline::setLambda (double lambda){
	if(offlineObject_->getConfig()->decomp_type_ == DBMatDecompEigen)
		offlineObject_->getConfig()->lambda_ = lambda;
	else if(offlineObject_->getConfig()->decomp_type_ == DBMatDecompChol)
		offlineObject_->getConfig()->lambda_ = lambda;
	else
		throw sgpp::base::application_exception("Lambda can not be changed in the online step for this decomposition type!");
}

void DBMatOnline::readOffline(DBMatOffline* o) {
  offlineObject_ = o;
}

DBMatOffline* DBMatOnline::getOffline() {
  return offlineObject_;
}
