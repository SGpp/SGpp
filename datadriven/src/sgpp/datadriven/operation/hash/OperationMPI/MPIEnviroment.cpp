// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "MPIOperationFactory.hpp"
#include "MPIEnviroment.hpp"

MPIEnviroment* MPIEnviroment::singleton_instance = NULL;

MPIEnviroment::MPIEnviroment(int argc, char *argv[],
                             bool verbose) : numTasks(0), rank(0), verbose(verbose) {
  MPI_Init(&argc, &argv);
  // Gets number of tasks/processes that this program is running on
  MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
  // Gets the rank (process/task number) that this program is running on
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}
MPIEnviroment::MPIEnviroment(void) {}
MPIEnviroment::MPIEnviroment(MPIEnviroment &cpy) {
}

void MPIEnviroment::slave_mainloop(void) {
  if (rank != 0) {
    MPI_Status stat;
    if (verbose) {
      std::cout << "Node " << rank << ": Started slave mainloop" << std::endl;
    }
    MPISlaveOperation* slave_op = NULL;
    do {
      int messagesize = 0;
      MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_INT, &messagesize);
      int *message = new int[messagesize];
      MPI_Recv(message, messagesize, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
               MPI_COMM_WORLD, &stat);
      if (messagesize != 1) {
        std::cout << "Warning! Received illegal message at node " << rank
                  << "! Going to ignore this..." << std::endl;
        delete [] message;
        continue;
      }
      if (message[0] == -1) {
        delete [] message;
        break;
      }
      if (message[0] == 1) {
        // Create operation here
        MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
        MPI_Get_count(&stat, MPI_CHAR, &messagesize);
        char *classname = new char[messagesize];
        MPI_Recv(classname, messagesize, MPI_CHAR, stat.MPI_SOURCE, stat.MPI_TAG,
                 MPI_COMM_WORLD, &stat);
        if (slave_op != NULL) {
          std::cout << "Warning! Overriding existing mpi operation with operation "
                    << classname << " in node " << rank << std::endl;
          delete slave_op;
        }
        slave_op = create_mpi_operation(classname);
        if (verbose) {
          std::cout << "Node " << rank << ": Created slave operation \""
                    << classname << "\"" << std::endl;
        }
      }
      if (message[0] == 2) {
        // run operation here
        if (verbose) {
          std::cout << "Node " << rank << ": Starting slave operation " << std::endl;
        }
        if (slave_op == NULL)
          throw std::logic_error("Trying to run an non existing slave operation!");
        slave_op->slave_code();
      }
      delete [] message;
    }while(true);
    MPI_Finalize();
    if (verbose) {
      std::cout << "Node " << rank << ": Exiting... " << std::endl;
    }
    std::exit(0);
  }
}

void MPIEnviroment::init(int argc, char *argv[], bool verbose) {
  if (singleton_instance == NULL) {
    singleton_instance = new MPIEnviroment(argc, argv, verbose);
    singleton_instance->slave_mainloop();
  } else {
    throw std::logic_error("Singleton class \"MPIEnviroment\" already initialized!");
  }
}

void MPIEnviroment::release(void) {
  if (singleton_instance != NULL) {
    int message[1];
    // Send the order to terminate all slaves processes!
    message[0] = -1;
    for (int i = 0; i < MPIEnviroment::get_node_count(); i++) {
      MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    delete singleton_instance;
  } else {
    throw std::logic_error("Singleton class \"MPIEnviroment\" not yet initialized!");
  }
}

int MPIEnviroment::get_node_rank(void) {
  if (singleton_instance != NULL)
    return singleton_instance->rank;
  else
    throw std::logic_error("Singleton class \"MPIEnviroment\" not yet initialized!");
}

int MPIEnviroment::get_node_count(void) {
  if (singleton_instance != NULL)
    return singleton_instance->numTasks;
  else
    throw std::logic_error("Singleton class \"MPIEnviroment\" not yet initialized!");
}

MPIEnviroment::~MPIEnviroment(void) {
}
