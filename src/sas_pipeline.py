#!/usr/bin/env python3
"""
SAS (Sufficient Alignment Support) Integrated Pipeline
Complete workflow for multi-technology sequencing data analysis

This script orchestrates the complete SAS pipeline integrating:
1. BE (Base Error) detection using short-read, HiFi, and ONT data
2. SE (Structural Error) detection using filtered ONT data

Features:
- Multi-technology data integration (Illumina, HiFi, ONT)
- Automatic reference index generation
- Optional BAM filtering skip
- Comprehensive error detection and classification
- Parallel processing support
- Detailed logging and reporting

License: MIT
"""

import os
import sys
import argparse
import subprocess
import logging
import shutil
from pathlib import Path
from datetime import datetime
import json
from concurrent.futures import ThreadPoolExecutor, as_completed


def get_pipeline_base_dir():
    """Get the base directory of the SAS Pipeline based on script location"""
    # Get the directory containing this script
    script_dir = Path(__file__).parent.resolve()
    
    # The script is in src/, so the base directory is the parent
    base_dir = script_dir.parent
    
    return base_dir


class SASPipeline:
    """Main pipeline class for SAS algorithm workflow"""
    
    def __init__(self, config):
        self.config = config
        self.output_dir = Path(config['output_dir'])
        
        # Create output directory first, before setting up logging
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Auto-detect pipeline base directory if not provided
        if config.get('scripts_dir') == 'scripts' or config.get('src_dir') == 'src':
            pipeline_base = get_pipeline_base_dir()
            self.scripts_dir = pipeline_base / 'scripts'
            self.src_dir = pipeline_base / 'src'
        else:
            self.scripts_dir = Path(config.get('scripts_dir', 'scripts'))
            self.src_dir = Path(config.get('src_dir', 'src'))
        
        # Setup logging after creating output directory
        self.logger = self._setup_logging()
        
        # Log the directories being used
        self.logger.info(f"Pipeline base directory: {get_pipeline_base_dir()}")
        self.logger.info(f"Scripts directory: {self.scripts_dir}")
        self.logger.info(f"Source directory: {self.src_dir}")
        
        # File paths - convert to absolute paths to avoid working directory issues
        self.short_read_bam = Path(config['short_read_bam']).resolve()
        self.hifi_bam = Path(config['hifi_bam']).resolve()
        self.ont_bam = Path(config['ont_bam']).resolve()
        self.reference_fa = Path(config['reference_fa']).resolve()
        self.reference_fai = Path(config.get('reference_fai')).resolve() if config.get('reference_fai') else None
        self.threads = config.get('threads', 8)
        
        # Log the resolved paths
        self.logger.info(f"Short-read BAM: {self.short_read_bam}")
        self.logger.info(f"HiFi BAM: {self.hifi_bam}")
        self.logger.info(f"ONT BAM: {self.ont_bam}")
        self.logger.info(f"Reference FASTA: {self.reference_fa}")
        if self.reference_fai:
            self.logger.info(f"Reference FAI: {self.reference_fai}")
        
        # Analysis options
        self.skip_filter = config.get('skip_filter', False)
        self.run_be_analysis = config.get('run_be_analysis', True)
        self.run_se_analysis = config.get('run_se_analysis', True)
        
        # SE analysis parameters
        self.se_window_size = config.get('se_window_size', 2000)
        self.se_min_indel_length = config.get('se_min_indel_length', 50)
        
        # Setup reference fai file
        self._setup_reference_fai()
        
        # Output file paths - adjust based on whether filtering is skipped
        if self.skip_filter:
            # Use input BAM files directly (already absolute paths)
            self.short_filtered_bam = self.short_read_bam
            self.hifi_filtered_bam = self.hifi_bam
            self.ont_filtered_bam = self.ont_bam
        else:
            # Use filtered output files
            self.short_filtered_bam = self.output_dir / "short_read_filtered.bam"
            self.hifi_filtered_bam = self.output_dir / "hifi_filtered.bam"
            self.ont_filtered_bam = self.output_dir / "ont_filtered.bam"
        
        # BE analysis output files
        self.be_output_dir = self.output_dir / "be_analysis"
        self.be_output_dir.mkdir(exist_ok=True)
        self.smallerror_bed = self.be_output_dir / "smallerror.bed"
        self.hifi_filtered_regions = self.be_output_dir / "hifi_filtered_regions.txt"
        self.ont_filtered_regions = self.be_output_dir / "ont_filtered_regions.txt"
        self.merged_results = self.be_output_dir / "merged_results.txt"
        self.final_be_results = self.be_output_dir / "final_be_results.txt"
        
        # SE analysis output files
        self.se_output_dir = self.output_dir / "se_analysis"
        self.se_output_dir.mkdir(exist_ok=True)
        self.se_results_prefix = self.se_output_dir / "structural_errors"
    
    def _setup_logging(self):
        """Setup logging configuration"""
        log_file = self.output_dir / "sas_pipeline.log"
        
        # Clear any existing handlers
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        
        handlers = [logging.StreamHandler(sys.stdout)]
        
        # Try to create file handler, but don't fail if it's not possible
        try:
            handlers.append(logging.FileHandler(log_file))
        except (OSError, IOError, PermissionError) as e:
            print(f"Warning: Could not create log file {log_file}: {e}", file=sys.stderr)
            print("Continuing with console logging only.", file=sys.stderr)
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=handlers
        )
        return logging.getLogger(__name__)
    
    def _setup_reference_fai(self):
        """Setup reference fai file, generate if not provided"""
        if self.reference_fai is None:
            # Generate fai file path using Path object methods
            self.reference_fai = Path(str(self.reference_fa) + ".fai")
            
        # Check if fai file exists
        if not self.reference_fai.exists():
            self.logger.info(f"Reference fai file not found: {self.reference_fai}")
            self.logger.info("Generating fai file using samtools faidx...")
            
            # Generate fai file using samtools faidx
            cmd = ["samtools", "faidx", str(self.reference_fa)]
            try:
                result = subprocess.run(
                    cmd,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                )
                self.logger.info(f"Successfully generated fai file: {self.reference_fai}")
                
            except subprocess.CalledProcessError as e:
                self.logger.error(f"Failed to generate fai file: {e}")
                self.logger.error(f"Command output: {e.stdout}")
                self.logger.error(f"Command error: {e.stderr}")
                raise
            except FileNotFoundError:
                self.logger.error("samtools not found. Please install samtools.")
                raise
        else:
            self.logger.info(f"Using existing fai file: {self.reference_fai}")
    
    def _run_command(self, cmd, description, cwd=None):
        """Run a command with proper error handling and logging"""
        self.logger.info(f"Starting: {description}")
        
        # Convert all command elements to strings for safe logging and execution
        cmd_str = [str(item) for item in cmd]
        self.logger.info(f"Command: {' '.join(cmd_str)}")
        
        try:
            result = subprocess.run(
                cmd_str,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                cwd=cwd  # 不设置默认的cwd，让子进程在当前工作目录运行
            )
            
            if result.stdout:
                self.logger.info(f"STDOUT: {result.stdout}")
            if result.stderr:
                self.logger.warning(f"STDERR: {result.stderr}")
                
            self.logger.info(f"Completed: {description}")
            return result
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed: {description}")
            self.logger.error(f"Return code: {e.returncode}")
            self.logger.error(f"STDOUT: {e.stdout}")
            self.logger.error(f"STDERR: {e.stderr}")
            raise
    
    def step1_filter_bam_files(self):
        """Step 1: Filter all three BAM files using appropriate modes"""
        self.logger.info("="*50)
        self.logger.info("STEP 1: Filtering BAM files")
        self.logger.info("="*50)
        
        if self.skip_filter:
            self.logger.info("Skipping BAM filtering - input files are already filtered")
            self.logger.info(f"Using input BAM files directly:")
            self.logger.info(f"  Short-read BAM: {self.short_read_bam}")
            self.logger.info(f"  HiFi BAM: {self.hifi_bam}")
            self.logger.info(f"  ONT BAM: {self.ont_bam}")
        else:
            # Filter short-read BAM (perfect matches, no clipping)
            cmd = [
                str(self.scripts_dir / "filter_bam.sh"),
                "short-read",
                "-i", str(self.short_read_bam),
                "-o", str(self.short_filtered_bam),
                "-t", str(self.threads)
            ]
            self._run_command(cmd, "Filtering short-read BAM")
            
            # Filter HiFi BAM (split alignment filtering)
            cmd = [
                str(self.scripts_dir / "filter_bam.sh"),
                "long-read",
                "-i", str(self.hifi_bam),
                "-o", str(self.hifi_filtered_bam),
                "-t", str(self.threads)
            ]
            self._run_command(cmd, "Filtering HiFi BAM")
            
            # Filter ONT BAM (split alignment filtering)
            cmd = [
                str(self.scripts_dir / "filter_bam.sh"),
                "long-read",
                "-i", str(self.ont_bam),
                "-o", str(self.ont_filtered_bam),
                "-t", str(self.threads)
            ]
            self._run_command(cmd, "Filtering ONT BAM")
            
            self.logger.info("Step 1 completed: All BAM files filtered successfully")
        
        self.logger.info("Step 1 completed: BAM file preparation finished")
    
    def step2_be_analysis(self):
        """Step 2: Run BE (Base Error) analysis"""
        if not self.run_be_analysis:
            self.logger.info("Skipping BE analysis as requested")
            return
            
        self.logger.info("="*50)
        self.logger.info("STEP 2: Base Error (BE) Analysis")
        self.logger.info("="*50)
        
        # Step 2a: Detect small errors from short-read data
        self.logger.info("Step 2a: Detecting small error regions from short-read data")
        output_prefix = str(self.be_output_dir / "smallerror_short")
        cmd = [
            str(self.scripts_dir / "smallerror_short.sh"),
            "-i", str(self.short_filtered_bam),
            "-o", output_prefix,
            "-t", str(self.threads)
        ]
        self._run_command(cmd, "Detecting small error regions from short-read data")
        
        # The output will be {output_prefix}.smallerror.bed
        expected_bed = Path(f"{output_prefix}.smallerror.bed")
        if expected_bed.exists():
            shutil.copy(expected_bed, self.smallerror_bed)
            self.logger.info(f"Small error BED file created: {self.smallerror_bed}")
        else:
            raise FileNotFoundError(f"Expected smallerror BED file not found: {expected_bed}")
        
        # Step 2b: Process with HiFi data (based on original smallerror.bed)
        self.logger.info("Step 2b: Processing candidate regions with HiFi data")
        cmd = [
            "python3",
            str(self.scripts_dir / "smallerror_long.py"),
            str(self.smallerror_bed),
            str(self.hifi_filtered_bam),
            str(self.reference_fa),
            str(self.hifi_filtered_regions),
            "--processes", str(self.threads)
        ]
        self._run_command(cmd, "Processing with HiFi data")
        
        # Step 2c: Process with ONT data (also based on original smallerror.bed)
        self.logger.info("Step 2c: Processing candidate regions with ONT data")
        cmd = [
            "python3",
            str(self.scripts_dir / "smallerror_long.py"),
            str(self.smallerror_bed),  # Use original smallerror.bed, not hifi_filtered_regions
            str(self.ont_filtered_bam),
            str(self.reference_fa),
            str(self.ont_filtered_regions),
            "--processes", str(self.threads)
        ]
        self._run_command(cmd, "Processing with ONT data")
        
        # Step 2d: Merge HiFi and ONT results (both processed from original smallerror.bed)
        self.logger.info("Step 2d: Merging HiFi and ONT analysis results")
        cmd = [
            "python3",
            str(self.scripts_dir / "merge_pileup.py"),
            str(self.smallerror_bed),
            str(self.hifi_filtered_regions),
            str(self.ont_filtered_regions),
            str(self.reference_fa),
            str(self.merged_results)
        ]
        self._run_command(cmd, "Merging HiFi and ONT pileup results")
        
        # Step 2e: Final BE filtering
        self.logger.info("Step 2e: Final BE filtering")
        cmd = [
            "python3",
            str(self.scripts_dir / "filter_be.py"),
            str(self.merged_results),
            str(self.reference_fai),
            "-o", str(self.final_be_results)
        ]
        self._run_command(cmd, "Final BE filtering")
        
        self.logger.info("BE analysis completed successfully")
    
    def step3_se_analysis(self):
        """Step 3: Run SE (Structural Error) analysis"""
        if not self.run_se_analysis:
            self.logger.info("Skipping SE analysis as requested")
            return
            
        self.logger.info("="*50)
        self.logger.info("STEP 3: Structural Error (SE) Analysis")
        self.logger.info("="*50)
        
        # Run integrated SE analysis using filtered ONT BAM
        cmd = [
            "python3",
            str(self.src_dir / "integrated_se_analysis.py"),
            str(self.ont_filtered_bam),
            str(self.reference_fa),
            "structural_errors",  # Only pass the prefix, not full path
            "--output-dir", str(self.se_output_dir),
            "--window-size", str(self.se_window_size),
            "--min-indel-length", str(self.se_min_indel_length),
            "--threads", str(self.threads)
        ]
        self._run_command(cmd, "Running SE analysis with filtered ONT data")
        
        # Update se_results_prefix to point to the actual output files
        self.se_results_prefix = self.se_output_dir / "structural_errors"
        
        self.logger.info("SE analysis completed successfully")
    
    def step4_generate_summary(self):
        """Step 4: Generate comprehensive summary report"""
        self.logger.info("="*50)
        self.logger.info("STEP 4: Generating Summary Report")
        self.logger.info("="*50)
        
        summary_file = self.output_dir / "sas_pipeline_summary.txt"
        
        # Count results
        def count_lines(file_path):
            try:
                with open(file_path, 'r') as f:
                    return sum(1 for _ in f)
            except:
                return 0
        
        # Get reference genome length from fai file
        def get_reference_length():
            try:
                total_length = 0
                with open(self.reference_fai, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            total_length += int(parts[1])
                return total_length
            except Exception as e:
                self.logger.warning(f"Could not read reference length from {self.reference_fai}: {e}")
                return None
        
        reference_length = get_reference_length()
        
        # BE analysis results
        be_results = {}
        total_be_errors = 0
        if self.run_be_analysis:
            be_results = {
                'smallerror_count': count_lines(self.smallerror_bed),
                'hifi_filtered_count': count_lines(self.hifi_filtered_regions) - 1,  # Exclude header
                'ont_filtered_count': count_lines(self.ont_filtered_regions) - 1,   # Exclude header
                'merged_count': count_lines(self.merged_results) - 1,               # Exclude header
                'final_be_count': count_lines(self.final_be_results)
            }
            total_be_errors = be_results['final_be_count']
        
        # SE analysis results
        se_results = {}
        total_se_errors = 0
        if self.run_se_analysis:
            se_anomalies_file = self.se_output_dir / "structural_errors.anomalous_regions.tsv"
            se_count = count_lines(se_anomalies_file) - 1 if se_anomalies_file.exists() else 0  # Exclude header
            se_results = {
                'anomalous_regions_count': se_count
            }
            total_se_errors = se_count
        
        # Calculate QV value
        total_errors = total_be_errors + total_se_errors
        qv_value = None
        qv_display = "N/A"
        
        if reference_length is not None:
            if total_errors == 0:
                qv_display = "∞ (Infinite - No errors detected)"
                self.logger.info("QV = ∞ (No errors detected)")
            else:
                import math
                qv_value = -10 * math.log10(total_errors / reference_length)
                qv_display = f"{qv_value:.2f}"
                self.logger.info(f"QV = {qv_value:.2f}")
        else:
            self.logger.warning("Cannot calculate QV: reference length unknown")
        
        # Store QV results as instance variables for final display
        self.final_qv_value = qv_value
        self.final_qv_display = qv_display
        self.final_total_errors = total_errors
        self.final_be_errors = total_be_errors
        self.final_se_errors = total_se_errors
        self.final_reference_length = reference_length

        with open(summary_file, 'w') as f:
            f.write("SAS Algorithm Pipeline Summary\n")
            f.write("="*50 + "\n")
            f.write(f"Execution Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Output Directory: {self.output_dir}\n")
            f.write(f"Threads Used: {self.threads}\n")
            f.write(f"Skip Filter: {self.skip_filter}\n")
            f.write(f"Run BE Analysis: {self.run_be_analysis}\n")
            f.write(f"Run SE Analysis: {self.run_se_analysis}\n\n")
            
            f.write("Input Files:\n")
            f.write(f"  Short-read BAM: {self.short_read_bam}\n")
            f.write(f"  HiFi BAM: {self.hifi_bam}\n")
            f.write(f"  ONT BAM: {self.ont_bam}\n")
            f.write(f"  Reference: {self.reference_fa}\n")
            f.write(f"  Reference FAI: {self.reference_fai}\n")
            if reference_length:
                f.write(f"  Reference Length: {reference_length:,} bp\n")
            f.write("\n")
            
            f.write("Output Files:\n")
            f.write(f"  Short-read filtered BAM: {self.short_filtered_bam}\n")
            f.write(f"  HiFi filtered BAM: {self.hifi_filtered_bam}\n")
            f.write(f"  ONT filtered BAM: {self.ont_filtered_bam}\n")
            
            if self.run_be_analysis:
                f.write(f"  BE Analysis Results: {self.be_output_dir}\n")
                f.write(f"  Final BE results: {self.final_be_results}\n")
            
            if self.run_se_analysis:
                f.write(f"  SE Analysis Results: {self.se_output_dir}\n")
                f.write(f"  SE anomalous regions: {self.se_output_dir}/structural_errors.anomalous_regions.tsv\n")
            
            f.write("\nResults Statistics:\n")
            
            if self.run_be_analysis:
                f.write("  BE (Base-level Errors) Analysis:\n")
                f.write(f"    Initial small error regions: {be_results['smallerror_count']}\n")
                f.write(f"    After HiFi filtering: {be_results['hifi_filtered_count']}\n")
                f.write(f"    After ONT filtering: {be_results['ont_filtered_count']}\n")
                f.write(f"    After merging: {be_results['merged_count']}\n")
                f.write(f"    Final BE count: {be_results['final_be_count']}\n")
            
            if self.run_se_analysis:
                f.write("  SE (Structural Errors) Analysis:\n")
                f.write(f"    Anomalous regions detected: {se_results['anomalous_regions_count']}\n")
            
            # Quality Value calculation
            f.write("\nQuality Assessment:\n")
            f.write(f"  Total BE errors: {total_be_errors}\n")
            f.write(f"  Total SE errors: {total_se_errors}\n")
            f.write(f"  Total errors: {total_errors}\n")
            if reference_length:
                f.write(f"  Assembly length: {reference_length:,} bp\n")
                f.write(f"  Error rate: {total_errors/reference_length:.2e}\n")
            f.write(f"  Quality Value (QV): {qv_display}\n")
            
            f.write(f"\nQV Formula: QV = -10 × log₁₀((SE + BE) / Assembly Length)\n")
            if total_errors == 0:
                f.write("Note: QV is infinite when no errors are detected.\n")
            
            f.write("\nPipeline Steps Completed:\n")
            if self.skip_filter:
                f.write("  ✓ Step 1: BAM files preparation (skipped filtering - using pre-filtered files)\n")
            else:
                f.write("  ✓ Step 1: BAM files filtered\n")
            
            if self.run_be_analysis:
                f.write("  ✓ Step 2: BE analysis completed\n")
            else:
                f.write("  - Step 2: BE analysis skipped\n")
            
            if self.run_se_analysis:
                f.write("  ✓ Step 3: SE analysis completed\n")
            else:
                f.write("  - Step 3: SE analysis skipped\n")
            
            f.write("  ✓ Step 4: Summary report generated\n")
        
        self.logger.info(f"Summary report saved: {summary_file}")
        
        # Log final QV to console for immediate feedback
        if total_errors == 0:
            self.logger.info("🎉 PERFECT ASSEMBLY: QV = ∞ (No errors detected)")
        elif qv_value is not None:
            self.logger.info(f"📊 FINAL QUALITY VALUE: QV = {qv_value:.2f}")
        
        # Also save configuration for reproducibility
        config_file = self.output_dir / "pipeline_config.json"
        with open(config_file, 'w') as f:
            json.dump(self.config, f, indent=2)
        
        self.logger.info(f"Configuration saved: {config_file}")
    
    def run_complete_pipeline(self):
        """Run the complete SAS pipeline"""
        start_time = datetime.now()
        
        self.logger.info("="*60)
        self.logger.info("SAS ALGORITHM INTEGRATED PIPELINE STARTED")
        self.logger.info("="*60)
        self.logger.info(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        try:
            self.step1_filter_bam_files()
            self.step2_be_analysis()
            self.step3_se_analysis()
            self.step4_generate_summary()
            
            end_time = datetime.now()
            duration = end_time - start_time
            
            self.logger.info("="*60)
            self.logger.info("SAS PIPELINE COMPLETED SUCCESSFULLY!")
            self.logger.info("="*60)
            self.logger.info(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
            self.logger.info(f"Total duration: {duration}")
            self.logger.info(f"All results saved in: {self.output_dir}")
            
            # Display final quality assessment
            if hasattr(self, 'final_qv_display'):
                self.logger.info("="*60)
                self.logger.info("FINAL QUALITY ASSESSMENT")
                self.logger.info("="*60)
                self.logger.info(f"📊 Total BE errors: {getattr(self, 'final_be_errors', 0)}")
                self.logger.info(f"📊 Total SE errors: {getattr(self, 'final_se_errors', 0)}")
                self.logger.info(f"📊 Total errors: {getattr(self, 'final_total_errors', 0)}")
                if hasattr(self, 'final_reference_length') and self.final_reference_length:
                    self.logger.info(f"📏 Assembly length: {self.final_reference_length:,} bp")
                    
                if self.final_total_errors == 0:
                    self.logger.info("🎉 ASSEMBLY QUALITY: QV = ∞ (No errors detected!)")
                elif self.final_qv_value is not None:
                    self.logger.info(f"📊 ASSEMBLY QUALITY: QV = {self.final_qv_value:.2f}")
                else:
                    self.logger.info(f"📊 ASSEMBLY QUALITY: QV = {self.final_qv_display}")
                self.logger.info("="*60)
            
            if self.run_be_analysis:
                self.logger.info(f"BE results: {self.final_be_results}")
            if self.run_se_analysis:
                self.logger.info(f"SE results: {self.se_output_dir}")
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            raise


def create_config_from_args(args):
    """Create configuration dictionary from command line arguments"""
    config = {
        'short_read_bam': args.short_read_bam,
        'hifi_bam': args.hifi_bam,
        'ont_bam': args.ont_bam,
        'reference_fa': args.reference_fa,
        'output_dir': args.output_dir,
        'threads': args.threads,
        'scripts_dir': args.scripts_dir,
        'src_dir': args.src_dir,
        'skip_filter': args.skip_filter,
        'run_be_analysis': args.run_be_analysis,
        'run_se_analysis': args.run_se_analysis,
        'se_window_size': args.se_window_size,
        'se_min_indel_length': args.se_min_indel_length
    }
    
    # Add reference_fai if provided
    if hasattr(args, 'reference_fai') and args.reference_fai:
        config['reference_fai'] = args.reference_fai
    
    return config


def validate_input_files(config):
    """Validate that all required input files exist"""
    required_files = [
        'short_read_bam',
        'hifi_bam', 
        'ont_bam',
        'reference_fa'
    ]
    
    for file_key in required_files:
        file_path = config[file_key]
        # Convert to Path object if it's a string
        if isinstance(file_path, str):
            file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Required input file not found: {file_path}")
    
    # Check reference_fai if provided
    if 'reference_fai' in config and config['reference_fai']:
        ref_fai_path = config['reference_fai']
        if isinstance(ref_fai_path, str):
            ref_fai_path = Path(ref_fai_path)
        if not ref_fai_path.exists():
            raise FileNotFoundError(f"Reference fai file not found: {ref_fai_path}")
    
    # Check scripts directory
    scripts_dir = Path(config['scripts_dir'])
    src_dir = Path(config['src_dir'])
    
    # Required scripts for BE analysis
    required_scripts = []
    if config.get('run_be_analysis', True):
        required_scripts.extend([
            'smallerror_short.sh',
            'smallerror_long.py',
            'merge_pileup.py',
            'filter_be.py'
        ])
        
        # Add filter.sh if not skipping filter
        if not config.get('skip_filter', False):
            required_scripts.append('filter_bam.sh')
    
    # Required scripts for SE analysis
    if config.get('run_se_analysis', True):
        required_src_files = ['integrated_se_analysis.py']
        for src_file in required_src_files:
            src_path = src_dir / src_file
            if not src_path.exists():
                raise FileNotFoundError(f"Required source file not found: {src_path}")
    
    for script in required_scripts:
        script_path = scripts_dir / script
        if not script_path.exists():
            raise FileNotFoundError(f"Required script not found: {script_path}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="SAS Algorithm Integrated Pipeline - Base Error and Structural Error Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline (BE + SE analysis)
  python3 sas_pipeline.py \\
    --short-read-bam illumina.bam \\
    --hifi-bam hifi.bam \\
    --ont-bam ont.bam \\
    --reference-fa reference.fa \\
    --output-dir results \\
    --threads 32

  # Run only BE analysis
  python3 sas_pipeline.py \\
    --short-read-bam illumina.bam \\
    --hifi-bam hifi.bam \\
    --ont-bam ont.bam \\
    --reference-fa reference.fa \\
    --output-dir results \\
    --no-se-analysis \\
    --threads 32

  # Run only SE analysis (with pre-filtered BAM)
  python3 sas_pipeline.py \\
    --short-read-bam illumina.bam \\
    --hifi-bam hifi.bam \\
    --ont-bam ont_filtered.bam \\
    --reference-fa reference.fa \\
    --output-dir results \\
    --skip-filter \\
    --no-be-analysis \\
    --threads 32

Output Structure:
  results/
    ├── be_analysis/          # Base error analysis results
    │   ├── final_be_results.txt
    │   └── ...
    ├── se_analysis/          # Structural error analysis results
    │   ├── structural_errors.anomalous_regions.tsv
    │   └── ...
    ├── sas_pipeline_summary.txt
    └── sas_pipeline.log
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--short-read-bam",
        required=True,
        help="Path to short-read BAM file (Illumina, Element, etc.)"
    )
    
    parser.add_argument(
        "--hifi-bam",
        required=True,
        help="Path to HiFi BAM file (PacBio HiFi)"
    )
    
    parser.add_argument(
        "--ont-bam",
        required=True,
        help="Path to ONT BAM file (Oxford Nanopore)"
    )
    
    parser.add_argument(
        "--reference-fa",
        required=True,
        help="Path to reference genome FASTA file"
    )
    
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for all results"
    )
    
    # Optional arguments
    parser.add_argument(
        "--reference-fai",
        help="Path to reference genome FASTA index file (.fai). If not provided, will be generated automatically"
    )
    
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads to use (default: 8)"
    )
    
    # Get auto-detected pipeline base directory
    pipeline_base = get_pipeline_base_dir()
    
    parser.add_argument(
        "--scripts-dir",
        default=str(pipeline_base / "scripts"),
        help=f"Path to scripts directory (default: {pipeline_base / 'scripts'})"
    )
    
    parser.add_argument(
        "--src-dir",
        default=str(pipeline_base / "src"),
        help=f"Path to source directory (default: {pipeline_base / 'src'})"
    )
    
    # Pipeline control options
    parser.add_argument(
        "--skip-filter",
        action="store_true",
        help="Skip BAM filtering step if input BAM files are already filtered"
    )
    
    parser.add_argument(
        "--no-be-analysis",
        action="store_true",
        help="Skip base error (BE) analysis"
    )
    
    parser.add_argument(
        "--no-se-analysis",
        action="store_true",
        help="Skip structural error (SE) analysis"
    )
    
    # SE analysis parameters
    parser.add_argument(
        "--se-window-size",
        type=int,
        default=2000,
        help="Window size for SE analysis (default: 2000)"
    )
    
    parser.add_argument(
        "--se-min-indel-length",
        type=int,
        default=50,
        help="Minimum INDEL length for SE analysis (default: 50)"
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="SAS Integrated Pipeline v2.0.0"
    )
    
    args = parser.parse_args()
    
    # Convert negative flags to positive
    args.run_be_analysis = not args.no_be_analysis
    args.run_se_analysis = not args.no_se_analysis
    
    # Validate that at least one analysis is enabled
    if not args.run_be_analysis and not args.run_se_analysis:
        print("Error: At least one analysis type (BE or SE) must be enabled", file=sys.stderr)
        sys.exit(1)
    
    # Create configuration
    config = create_config_from_args(args)
    
    # Validate input files
    try:
        validate_input_files(config)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Run pipeline
    try:
        pipeline = SASPipeline(config)
        pipeline.run_complete_pipeline()
        
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Pipeline failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 