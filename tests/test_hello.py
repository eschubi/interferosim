"""Tests for the hello-world example."""

from interferosim.hello import hello_world, main


def test_hello_world_default():
    assert hello_world() == "Hello, world!"


def test_hello_world_custom_name():
    assert hello_world("Interferosim") == "Hello, Interferosim!"


def test_main_prints_greeting(capsys):
    exit_code = main(["Interferosim"])
    captured = capsys.readouterr()
    assert exit_code == 0
    assert captured.out.strip() == "Hello, Interferosim!"
