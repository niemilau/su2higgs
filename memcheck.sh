#!/bin/bash

#best to test this with serial code. MPI seems to lead to tons of errors..

valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out ./su2 config

