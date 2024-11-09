# -*- coding: utf-8 -*-
# @Time    : 2024/10/8 21:43
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : biotool.py

import os
import subprocess


class BioCommandWrapper:
    def __init__(self, command_name, version_command=None):
        if not version_command:
            version_command = ['--version']
        self.command_name = command_name
        self.version_command = version_command
        self.check_command()

    def check_command(self):
        try:
            result = subprocess.run([self.command_name] + self.version_command, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            if result.returncode != 0:
                raise FileNotFoundError(f"Can't find the {self.command_name} executable.")
        except FileNotFoundError:
            raise FileNotFoundError(f"Can't find the {self.command_name} executable.")
        except Exception as e:
            raise Exception(f"Error occurred when checking the command: {self.command_name}\n{e}")

    def run_command(self, params: list):
        cmd = [self.command_name] + params
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode != 0:
                raise Exception(f"Command failed with error: {result.stderr.decode('utf-8')}")
            return result.stdout.decode('utf-8')
        except FileNotFoundError:
            raise FileNotFoundError(f"Can't find the executable for the command: {' '.join(cmd)}")
        except Exception as e:
            raise Exception(f"Error occurred when running the command: {' '.join(cmd)}\n{e}")

    def get_version(self):
        try:
            result = subprocess.run([self.command_name] + self.version_command, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            if result.returncode != 0:
                raise Exception(f"Failed to get version for {self.command_name}")
            return result.stdout.decode('utf-8').strip()
        except Exception as e:
            raise Exception(f"Error occurred when getting the version: {self.command_name}\n{e}")


class GemmaPipeline:
    def __init__(self, gemma_path='gemma', plink_path='plink'):
        self.gemma = BioCommandWrapper(gemma_path, ['-license'])
        self.plink = BioCommandWrapper(plink_path, ['--version'])

        self.gemma.check_command()
        self.plink.check_command()

    def run(self, vcf_file=None, phe_file=None, out_name="gwas", out_dir=".", phe_num=1):
        """
        Run the GEMMA analysis for all phenotype
        :param vcf_file: a VCF file containing the genotype data
        :param phe_file: a phenotype file containing the phenotype data, can be multiple phenotype
        :param out_name: the output file name, default is "gwas"
        :param out_dir: the output directory, default is "."
        :param phe_num: the number of phenotype in the phenotype file, default is 1
        :return:
        """
        # Ensure out_dir ends with a slash
        if not out_dir.endswith('/'):
            out_dir += '/'

        # Check whether the input files exist
        if not os.path.exists(vcf_file):
            raise FileNotFoundError(f"The VCF file {vcf_file} does not exist.")
        if not os.path.exists(phe_file):
            raise FileNotFoundError(f"The phenotype file {phe_file} does not exist.")

        # Run the GEMMA analysis for all phenotype
        for i in range(1, phe_num + 1):
            # Full path for output files
            full_out_name = f"{out_dir}{out_name}_{i}"

            # Convert vcf file into plink format
            self.plink.run_command([
                '--vcf', vcf_file, '--pheno', phe_file, '--mpheno', str(i),
                '--out', full_out_name, '--missing-phenotype', '-9',
                '--vcf-half-call', 'm', '--allow-extra-chr', '--allow-no-sex', '--make-bed'
            ])
            # Estimate Relatedness Matrix from Genotypes
            self.gemma.run_command(
                ['-bfile', full_out_name, '-gk', '2', '-o', f"{out_name}_{i}", '-outdir', out_dir]
            )
            # Perform Association Tests with Univariate Linear Mixed Models
            self.gemma.run_command(
                ['-bfile', full_out_name, '-k', f'{full_out_name}.sXX.txt', '-lmm', '4', '-o',
                 f"{out_name}_{i}", '-outdir', out_dir]
            )
