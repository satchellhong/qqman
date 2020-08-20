import setuptools

with open("README.md", 'r') as fh:
    long_description = fh.read()
	
setuptools.setup(
    name="qqman",
    version="1.0.2",
    license='MIT',
    author="chol hong",
    author_email="shulkhorn@gmail.com",
    description="Draws Manhattan plot and QQ plot using plink assoc output.",
    long_description=long_description,
    long_description_content_type="text/markdown",
	url="https://github.com/satchellhong/qqman",
    packages=setuptools.find_packages(),
    classifiers=[
        # 패키지에 대한 태그
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
)