from setuptools import find_packages, setup

setup(
    name="src",
    version="0.0.1",
    author="Simon",
    author_email="author@example.com",
    description="Python repository for simulating particles moving around under gravity",
    url="url-to-github-page",
    packages=find_packages(),
    test_suite="src.tests.test_all.suite",
)
