from setuptools import setup

setup(
    name="raccoon",
    version="1.0.0",
    description="A Python program to generate 3D structures of linear polymers",
    author="Moritz L. Obenauer",
    url="https://github.com/moritzobenauer/ProjectRaccoon",
    packages=["raccoon"],
    package_data={"raccoon": ["src/data/*.json"]},
    entry_points={"console_scripts": ["raccoon=raccoon.__main__:main"]},
    install_requires=[
        "biopandas==0.4.1",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
