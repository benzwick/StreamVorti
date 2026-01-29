/**
 * @file test_main.cpp
 * @brief Main entry point for StreamVorti unit tests
 *
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2020-2025 Benjamin F. Zwick
 */

#include <gtest/gtest.h>

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
