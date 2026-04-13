/**
 * @file timer_output.hpp
 * @brief Section-based timer with formatted summary output.
 *
 * Inspired by deal.II's TimerOutput class. Provides named sections that
 * accumulate wall/CPU time across multiple enter/leave calls, with RAII
 * scoped timers and formatted summary tables.
 *
 * @see https://dealii.org/current/doxygen/deal.II/classTimerOutput.html
 */

#ifndef STREAMVORTI_TIMER_OUTPUT_HPP
#define STREAMVORTI_TIMER_OUTPUT_HPP

#include <mfem.hpp>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace StreamVorti {

/**
 * @brief Section-based timer with formatted summary output.
 *
 * Usage:
 * @code
 *   TimerOutput timer(std::cout, TimerOutput::summary,
 *                     TimerOutput::wall_times);
 *
 *   {
 *       TimerOutput::Scope t(timer, "Assembly");
 *       assemble_system();
 *   }
 *
 *   {
 *       TimerOutput::Scope t(timer, "Solve");
 *       solve_system();
 *   }
 *   // Summary printed automatically at destruction
 * @endcode
 *
 * MPI usage (with barrier synchronization):
 * @code
 *   TimerOutput timer(MPI_COMM_WORLD, std::cout,
 *                     TimerOutput::summary, TimerOutput::wall_times);
 * @endcode
 */
class TimerOutput
{
public:
    // ================================================================
    // Enumerations
    // ================================================================

    /// When to produce output.
    enum OutputFrequency
    {
        every_call,              ///< Print after every leave_subsection()
        summary,                 ///< Print summary at destruction
        every_call_and_summary,  ///< Both
        never                    ///< No automatic output
    };

    /// What data to return from get_summary_data().
    enum OutputData
    {
        total_cpu_time,
        total_wall_time,
        n_calls
    };

    /// What kind of timing to display.
    enum OutputType
    {
        cpu_times,                  ///< CPU times only
        wall_times,                 ///< Wall clock times only
        cpu_and_wall_times,         ///< Both in separate tables
        cpu_and_wall_times_grouped  ///< Both in a single table
    };

    // ================================================================
    // Nested types
    // ================================================================

    /**
     * @brief RAII scoped timer.
     *
     * Enters a section on construction, leaves it on destruction.
     * Exception-safe: guarantees leave even during stack unwinding.
     */
    class Scope
    {
    public:
        Scope(TimerOutput &timer, const std::string &section_name)
            : timer_(timer), section_name_(section_name), active_(true)
        {
            timer_.enter_subsection(section_name_);
        }

        ~Scope()
        {
            try { stop(); } catch (...) {}
        }

        /// Manually stop the timer before the scope ends.
        void stop()
        {
            if (active_)
            {
                timer_.leave_subsection(section_name_);
                active_ = false;
            }
        }

        Scope(const Scope &) = delete;
        Scope &operator=(const Scope &) = delete;

    private:
        TimerOutput &timer_;
        const std::string section_name_;
        bool active_;
    };

    // ================================================================
    // Constructors
    // ================================================================

    /// Serial constructor with std::ostream.
    TimerOutput(std::ostream &stream,
                const OutputFrequency output_frequency,
                const OutputType output_type)
        : output_frequency_(output_frequency)
        , output_type_(output_type)
        , stream_(&stream)
        , output_is_enabled_(true)
        , is_rank_zero_(true)
#ifdef MFEM_USE_MPI
        , mpi_communicator_(MPI_COMM_SELF)
        , has_mpi_(false)
#endif
    {
        timer_all_.Clear();
        timer_all_.Start();
    }

#ifdef MFEM_USE_MPI
    /// MPI constructor. Inserts barriers for synchronized timing.
    /// Only rank 0 produces output.
    TimerOutput(MPI_Comm mpi_comm,
                std::ostream &stream,
                const OutputFrequency output_frequency,
                const OutputType output_type)
        : output_frequency_(output_frequency)
        , output_type_(output_type)
        , stream_(&stream)
        , output_is_enabled_(true)
        , mpi_communicator_(mpi_comm)
        , has_mpi_(mpi_comm != MPI_COMM_SELF)
    {
        int rank = 0;
        MPI_Comm_rank(mpi_comm, &rank);
        is_rank_zero_ = (rank == 0);

        timer_all_.Clear();
        timer_all_.Start();
    }
#endif

    // ================================================================
    // Destructor
    // ================================================================

    ~TimerOutput()
    {
        auto do_exit = [this]() {
            try
            {
                while (!active_sections_.empty())
                    leave_subsection();

                if ((output_frequency_ == summary ||
                     output_frequency_ == every_call_and_summary) &&
                    output_is_enabled_)
                {
                    print_summary();
                }
            }
            catch (...) {}
        };

#ifdef MFEM_USE_MPI
        if (std::uncaught_exceptions() > 0 && has_mpi_)
        {
            if (is_rank_zero_)
            {
                std::cerr
                    << "---------------------------------------------------------\n"
                    << "TimerOutput: Skipping summary output during exception\n"
                    << "unwinding to avoid MPI deadlock.\n"
                    << "---------------------------------------------------------"
                    << std::endl;
            }
        }
        else
        {
            do_exit();
        }
#else
        do_exit();
#endif
    }

    // Non-copyable
    TimerOutput(const TimerOutput &) = delete;
    TimerOutput &operator=(const TimerOutput &) = delete;

    // ================================================================
    // Section entry/exit
    // ================================================================

    /// Enter a named section. Creates it if it doesn't exist.
    void enter_subsection(const std::string &section_name)
    {
        std::lock_guard<std::mutex> lock(mutex_);

        assert(!section_name.empty() && "Section name must not be empty");
        assert(std::find(active_sections_.begin(), active_sections_.end(),
                         section_name) == active_sections_.end() &&
               "Cannot enter an already active section");

#ifdef MFEM_USE_MPI
        if (has_mpi_)
            MPI_Barrier(mpi_communicator_);
#endif

        auto it = sections_.find(section_name);
        if (it == sections_.end())
        {
            // New section — record insertion order and parent
            section_order_.push_back(section_name);
            Section s;
            if (!active_sections_.empty())
            {
                s.parent = active_sections_.back();
                s.depth = sections_.at(s.parent).depth + 1;
            }
            sections_.emplace(section_name, std::move(s));
        }

        sections_.at(section_name).start();
        active_sections_.push_back(section_name);
    }

    /// Leave a named section. If empty, leaves the most recently entered.
    void leave_subsection(const std::string &section_name = "")
    {
        assert(!active_sections_.empty() &&
               "Cannot leave: no section is active");

        std::lock_guard<std::mutex> lock(mutex_);

#ifdef MFEM_USE_MPI
        if (has_mpi_)
            MPI_Barrier(mpi_communicator_);
#endif

        const std::string &actual_name =
            section_name.empty() ? active_sections_.back() : section_name;

        assert(sections_.count(actual_name) &&
               "Cannot leave a section that was never entered");
        assert(sections_[actual_name].running &&
               "Cannot leave a section that is not running");

        sections_[actual_name].stop();

        // Print per-call output
        if ((output_frequency_ == every_call ||
             output_frequency_ == every_call_and_summary) &&
            output_is_enabled_ && is_rank_zero_ &&
            std::uncaught_exceptions() == 0)
        {
            const Section &s = sections_[actual_name];
            *stream_ << actual_name;

            if (output_type_ == cpu_times)
                *stream_ << ", CPU time: " << s.last_cpu_time << "s";
            else if (output_type_ == wall_times)
                *stream_ << ", wall time: " << s.last_wall_time << "s.";
            else
                *stream_ << ", CPU/wall time: " << s.last_cpu_time
                         << "s / " << s.last_wall_time << "s.";

            *stream_ << std::endl;
        }

        // Remove from active list
        auto it = std::find(active_sections_.begin(),
                            active_sections_.end(), actual_name);
        assert(it != active_sections_.end());
        active_sections_.erase(it);
    }

    // ================================================================
    // Output
    // ================================================================

    /// Print a formatted summary table of all sections.
    void print_summary() const
    {
        if (!is_rank_zero_ || !output_is_enabled_)
            return;

        assert(active_sections_.empty() &&
               "Cannot print summary while inside a timed section");

        // Save stream state
        std::ios old_state(nullptr);
        old_state.copyfmt(*stream_);

        // Compute column width from longest section name (including indentation)
        unsigned int max_width = 32;
        for (const auto &name : section_order_)
        {
            unsigned int depth = sections_.at(name).depth;
            max_width = std::max(max_width,
                                 static_cast<unsigned int>(name.size() + depth * 2) + 1);
        }

        const std::string extra_dash(max_width - 32, '-');
        const std::string extra_space(max_width - 32, ' ');

        const double total_wall = timer_all_.RealTime();
        const double total_cpu = timer_all_.UserTime();

        if (output_type_ != cpu_and_wall_times_grouped)
        {
            // --- CPU times table ---
            if (output_type_ != wall_times)
            {
                double effective_total = total_cpu;
                double sum_cpu = 0;
                for (const auto &kv : sections_)
                    sum_cpu += kv.second.total_cpu_time;
                if (sum_cpu > effective_total)
                    effective_total = sum_cpu;

                print_table("CPU time", "Total CPU time elapsed since start",
                            effective_total, extra_dash, extra_space,
                            max_width, true);

                if (sum_cpu - total_cpu > 0.0)
                    print_overhead_note(sum_cpu - total_cpu);
            }

            // --- Wall times table ---
            if (output_type_ != cpu_times)
            {
                print_table("wall time",
                            "Total wallclock time elapsed since start",
                            total_wall, extra_dash, extra_space,
                            max_width, false);
            }
        }
        else
        {
            // --- Grouped CPU + wall table ---
            double effective_cpu = total_cpu;
            double sum_cpu = 0;
            for (const auto &kv : sections_)
                sum_cpu += kv.second.total_cpu_time;
            if (sum_cpu > effective_cpu)
                effective_cpu = sum_cpu;

            print_grouped_table(effective_cpu, total_wall,
                                extra_dash, extra_space, max_width);

            if (sum_cpu - total_cpu > 0.0)
                print_overhead_note(sum_cpu - total_cpu);
        }

        // Restore stream state
        stream_->copyfmt(old_state);
    }

#ifdef MFEM_USE_MPI
    /**
     * @brief Print wall time statistics across MPI ranks.
     *
     * Shows min, avg, max (and optionally quantiles) for each section.
     * Useful when the TimerOutput was constructed without an MPI communicator
     * (no barriers), to reveal load imbalance.
     *
     * @param mpi_comm  MPI communicator for gathering statistics.
     * @param quantile  If > 0, also print this quantile and (1-quantile).
     *                  E.g., 0.1 prints the 10th and 90th percentiles.
     */
    void print_wall_time_statistics(MPI_Comm mpi_comm,
                                    double quantile = 0.0) const
    {
        if (!output_is_enabled_)
            return;

        assert(active_sections_.empty() &&
               "Cannot print statistics while inside a timed section");

        int rank, nprocs;
        MPI_Comm_rank(mpi_comm, &rank);
        MPI_Comm_size(mpi_comm, &nprocs);

        // Save stream state
        std::ios old_state(nullptr);
        if (is_rank_zero_)
            old_state.copyfmt(*stream_);

        // Compute column width (including indentation)
        unsigned int max_width = 32;
        for (const auto &name : section_order_)
        {
            auto it = sections_.find(name);
            unsigned int depth = (it != sections_.end()) ? it->second.depth : 0;
            max_width = std::max(max_width,
                                 static_cast<unsigned int>(name.size() + depth * 2) + 1);
        }

        const std::string extra_dash(max_width - 32, '-');
        const std::string extra_space(max_width - 32, ' ');

        // Gather total wall time
        double local_total = timer_all_.RealTime();
        double global_max_total;
        MPI_Reduce(&local_total, &global_max_total, 1, MPI_DOUBLE,
                   MPI_MAX, 0, mpi_comm);

        bool print_quantiles = (quantile > 0.0 && quantile <= 0.5);

        if (is_rank_zero_)
        {
            *stream_ << "\n\n"
                      << "+---------------------------------------------"
                      << extra_dash;
            if (print_quantiles)
                *stream_ << "+------------+------------"
                         << "+------------+------------"
                         << "+------------+\n";
            else
                *stream_ << "+------------+------------"
                         << "+------------+\n";

            *stream_ << "| Total wallclock time elapsed since start    "
                      << extra_space << "|"
                      << std::setw(10) << std::setprecision(3) << std::right
                      << global_max_total << "s |"
                      << std::string(print_quantiles ? 49 : 25, ' ')
                      << "|\n";

            *stream_ << "|                                             "
                      << extra_space << "|"
                      << std::string(print_quantiles ? 61 : 37, ' ')
                      << "|\n";

            *stream_ << "| Section                         " << extra_space
                      << "| no. calls |"
                      << "    min (s) |    avg (s) |";
            if (print_quantiles)
                *stream_ << " " << std::setw(5) << std::setprecision(1)
                         << std::fixed << quantile * 100 << "% (s) |"
                         << " " << std::setw(5) << (1 - quantile) * 100
                         << "% (s) |";
            *stream_ << "    max (s) |\n";

            *stream_ << "+---------------------------------" << extra_dash
                      << "+-----------+";
            *stream_ << "------------+------------+";
            if (print_quantiles)
                *stream_ << "------------+------------+";
            *stream_ << "------------+";
        }

        // For each section, gather statistics
        for (const auto &name : section_order_)
        {
            auto it = sections_.find(name);
            double local_time = (it != sections_.end()) ?
                                it->second.total_wall_time : 0.0;
            unsigned int local_laps = (it != sections_.end()) ?
                                      it->second.n_laps : 0;

            // Gather all times to rank 0
            std::vector<double> all_times(nprocs, 0.0);
            MPI_Gather(&local_time, 1, MPI_DOUBLE,
                       all_times.data(), 1, MPI_DOUBLE, 0, mpi_comm);

            unsigned int global_laps = 0;
            MPI_Reduce(&local_laps, &global_laps, 1, MPI_UNSIGNED,
                       MPI_MAX, 0, mpi_comm);

            if (is_rank_zero_)
            {
                std::sort(all_times.begin(), all_times.end());
                double min_time = all_times.front();
                double max_time = all_times.back();
                double avg_time = 0;
                for (double t : all_times) avg_time += t;
                avg_time /= nprocs;

                unsigned int depth = (it != sections_.end()) ?
                                     it->second.depth : 0;
                std::string name_out = std::string(depth * 2, ' ') + name;
                name_out.resize(max_width, ' ');

                *stream_ << "\n| " << name_out << "| "
                         << std::setw(9) << global_laps << " |"
                         << std::setw(10) << std::setprecision(3)
                         << std::fixed << min_time << "s |"
                         << std::setw(10) << avg_time << "s |";

                if (print_quantiles)
                {
                    int lo_idx = static_cast<int>(
                        std::round(quantile * (nprocs - 1)));
                    int hi_idx = static_cast<int>(
                        std::round((1.0 - quantile) * (nprocs - 1)));
                    *stream_ << std::setw(10) << all_times[lo_idx] << "s |"
                             << std::setw(10) << all_times[hi_idx] << "s |";
                }

                *stream_ << std::setw(10) << max_time << "s |";
            }
        }

        if (is_rank_zero_)
        {
            *stream_ << "\n+---------------------------------" << extra_dash
                      << "+-----------+";
            *stream_ << "------------+------------+";
            if (print_quantiles)
                *stream_ << "------------+------------+";
            *stream_ << "------------+\n" << std::endl;

            stream_->copyfmt(old_state);
        }
    }
#endif

    /// Return a map of section_name -> value for the requested data type.
    std::map<std::string, double>
    get_summary_data(const OutputData kind) const
    {
        std::map<std::string, double> result;
        for (const auto &kv : sections_)
        {
            switch (kind)
            {
            case total_cpu_time:
                result[kv.first] = kv.second.total_cpu_time;
                break;
            case total_wall_time:
                result[kv.first] = kv.second.total_wall_time;
                break;
            case n_calls:
                result[kv.first] = static_cast<double>(kv.second.n_laps);
                break;
            }
        }
        return result;
    }

    /// Disable all output (does not affect timing).
    void disable_output() { output_is_enabled_ = false; }

    /// Re-enable output.
    void enable_output() { output_is_enabled_ = true; }

    /// Clear all sections and restart the overall timer.
    void reset()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        sections_.clear();
        section_order_.clear();
        active_sections_.clear();
        timer_all_.Clear();
        timer_all_.Start();
    }

private:
    // ================================================================
    // Internal section data
    // ================================================================

    struct Section
    {
        double total_wall_time = 0.0;
        double total_cpu_time = 0.0;
        double last_wall_time = 0.0;
        double last_cpu_time = 0.0;
        unsigned int n_laps = 0;
        bool running = false;
        std::unique_ptr<mfem::StopWatch> timer;  ///< Used for current lap
        std::string parent;                       ///< Parent section name ("" = top-level)
        unsigned int depth = 0;                   ///< Nesting depth (0 = top-level)

        Section() : timer(std::make_unique<mfem::StopWatch>()) {}

        void start()
        {
            timer->Clear();
            timer->Start();
            running = true;
        }

        void stop()
        {
            timer->Stop();
            last_wall_time = timer->RealTime();
            last_cpu_time = timer->UserTime();
            total_wall_time += last_wall_time;
            total_cpu_time += last_cpu_time;
            ++n_laps;
            running = false;
        }
    };

    // ================================================================
    // Print helpers
    // ================================================================

    /// Get the reference time for computing a section's percentage.
    /// Top-level sections use total_time; subsections use their parent's time.
    double get_reference_time(const Section &s, double total_time,
                              bool use_cpu) const
    {
        if (s.parent.empty())
            return total_time;
        const Section &parent = sections_.at(s.parent);
        return use_cpu ? parent.total_cpu_time : parent.total_wall_time;
    }

    /// Helper: format a percentage into a stream field.
    void print_percentage(double section_time, double ref_time) const
    {
        *stream_ << std::setw(10);
        if (ref_time != 0)
        {
            double fraction = section_time / ref_time;
            if (fraction > 0.001)
                *stream_ << std::setprecision(2) << fraction * 100;
            else
                *stream_ << 0.0;
            *stream_ << "% |";
        }
        else
            *stream_ << 0.0 << "% |";
    }

    /// Check if any section has a parent (i.e., subsections exist).
    bool has_subsections() const
    {
        for (const auto &name : section_order_)
            if (sections_.at(name).depth > 0)
                return true;
        return false;
    }

    /// Print a single-metric table (CPU or wall).
    void print_table(const std::string &time_label,
                     const std::string &header_text,
                     double total_time,
                     const std::string &extra_dash,
                     const std::string &extra_space,
                     unsigned int max_width,
                     bool use_cpu) const
    {
        *stream_ << "\n\n"
                  << "+---------------------------------------------"
                  << extra_dash << "+------------+------------+\n";

        // Pad header_text to fill the column
        std::string header_padded = header_text;
        header_padded.resize(44 + extra_space.size(), ' ');
        *stream_ << "| " << header_padded << "|"
                  << std::setw(10) << std::setprecision(3) << std::right
                  << total_time << "s |            |\n";

        *stream_ << "|                                             "
                  << extra_space << "|            |            |\n";

        // Column headers
        std::string section_header = "Section";
        section_header.resize(max_width, ' ');
        *stream_ << "| " << section_header << "| no. calls |"
                  << "  " << std::setw(9) << std::left << time_label
                  << std::right
                  << " | % of total |\n";

        *stream_ << "+---------------------------------" << extra_dash
                  << "+-----------+------------+------------+";

        for (const auto &name : section_order_)
        {
            const Section &s = sections_.at(name);
            double section_time = use_cpu ? s.total_cpu_time
                                          : s.total_wall_time;

            // Indent name by depth (2 spaces per level)
            std::string name_out = std::string(s.depth * 2, ' ') + name;
            name_out.resize(max_width, ' ');

            *stream_ << "\n| " << name_out << "| "
                      << std::setw(9) << s.n_laps << " |"
                      << std::setw(10) << std::setprecision(3) << std::fixed
                      << section_time << "s |";

            print_percentage(section_time, total_time);
        }

        *stream_ << "\n+---------------------------------" << extra_dash
                  << "+-----------+------------+------------+\n";

        if (has_subsections())
            *stream_ << "Indented sections show % of total;"
                      << " subsections of the same parent sum to ~parent %.\n";

        *stream_ << std::endl;
    }

    /// Print a grouped CPU + wall table.
    void print_grouped_table(double total_cpu, double total_wall,
                             const std::string &extra_dash,
                             const std::string &extra_space,
                             unsigned int max_width) const
    {
        // Top border
        *stream_ << "\n\n+---------------------------------------------"
                  << extra_dash << "+"
                  << "------------+------------+"
                  << "------------+------------+\n";

        // Header row
        std::string header = "Total CPU/wall time elapsed since start";
        header.resize(44 + extra_space.size(), ' ');
        *stream_ << "| " << header << "|"
                  << std::setw(10) << std::setprecision(3) << std::right
                  << total_cpu << "s |            |"
                  << std::setw(10) << total_wall << "s |            |\n";

        // Blank row
        *stream_ << "|                                             "
                  << extra_space << "|"
                  << "            |            |"
                  << "            |            |\n";

        // Column headers
        std::string section_header = "Section";
        section_header.resize(max_width, ' ');
        *stream_ << "| " << section_header << "| no. calls |"
                  << "  CPU time  | % of total |"
                  << "  wall time | % of total |\n";

        // Separator
        *stream_ << "+---------------------------------" << extra_dash
                  << "+-----------+"
                  << "------------+------------+"
                  << "------------+------------+" << std::endl;

        for (const auto &name : section_order_)
        {
            const Section &s = sections_.at(name);
            std::string name_out = std::string(s.depth * 2, ' ') + name;
            name_out.resize(max_width, ' ');

            *stream_ << "| " << name_out << "| "
                      << std::setw(9) << s.n_laps << " |";

            // CPU time + % of total
            *stream_ << std::setw(10) << std::setprecision(3) << std::fixed
                      << s.total_cpu_time << "s |";
            print_percentage(s.total_cpu_time, total_cpu);

            // Wall time + % of total
            *stream_ << std::setw(10) << std::setprecision(3) << std::fixed
                      << s.total_wall_time << "s |";
            print_percentage(s.total_wall_time, total_wall);

            *stream_ << "\n";
        }

        // Bottom border
        *stream_ << "+---------------------------------" << extra_dash
                  << "+-----------+"
                  << "------------+------------+"
                  << "------------+------------+\n";

        if (has_subsections())
            *stream_ << "Indented sections show % of total;"
                      << " subsections of the same parent sum to ~parent %.\n";

        *stream_ << std::endl;
    }

    void print_overhead_note(double overhead) const
    {
        *stream_ << "\nNote: The sum of section times is " << overhead
                  << " seconds larger than the total time.\n"
                  << "(Timer overhead, or sections may have run concurrently.)\n"
                  << std::endl;
    }

    // ================================================================
    // Member data
    // ================================================================

    OutputFrequency output_frequency_;
    OutputType output_type_;
    mutable mfem::StopWatch timer_all_;             ///< Overall elapsed timer
    std::map<std::string, Section> sections_;       ///< Section data by name
    std::vector<std::string> section_order_;         ///< Insertion order
    std::ostream *stream_;
    bool output_is_enabled_;
    bool is_rank_zero_;
    std::list<std::string> active_sections_;        ///< Stack of active sections
    std::mutex mutex_;

#ifdef MFEM_USE_MPI
    MPI_Comm mpi_communicator_;
    bool has_mpi_;                                  ///< True if comm != SELF
#endif
};

} // namespace StreamVorti

#endif // STREAMVORTI_TIMER_OUTPUT_HPP
